#include <checkRepeats.hpp>

namespace {

bool copyNumberInsideThreshold(float cn, float med, float th)
{
    if (cn <= med + th && cn >= med - th) {
        return cn > 0.1f;
    }
    return false;
}

SEGMENTINFO createSegments(
    const CONTIGLIST& contigs,
    uint32_t k_len,
    uint32_t count_threshold,
    float coverage_depth)
{
    SEGMENTINFO valSegment;
    if (contigs.empty()) {
        return valSegment;
    }

    std::string currentContig = contigs[0].contigName;
    uint32_t segStartBase = contigs[0].start;
    uint32_t segEndBase = contigs[0].end;
    uint32_t kCount = 1;
    const uint32_t gap_len = k_len;
    float wMedian = static_cast<float>(contigs[0].numRepeats) / coverage_depth;

    for (size_t i = 1; i < contigs.size(); i++)
    {
        if (copyNumberInsideThreshold(static_cast<float>(contigs[i].numRepeats) / coverage_depth, wMedian / kCount, static_cast<float>(count_threshold)) &&
            (contigs[i].contigName == currentContig) &&
            (contigs[i].end >= segEndBase) &&
            (contigs[i].end - segEndBase <= gap_len))
        {
            segEndBase = contigs[i].end;
            wMedian += static_cast<float>(contigs[i].numRepeats) / coverage_depth;
            kCount += 1;
        }
        else
        {
            if (contigs[i].start < segEndBase && contigs[i].contigName == currentContig)
            {
                if (segEndBase > segStartBase) {
                    valSegment.emplace_back(currentContig, segStartBase, segEndBase, wMedian / kCount);
                    segStartBase = segEndBase;
                }
            }
            else
            {
                if (segEndBase > segStartBase) {
                    valSegment.emplace_back(currentContig, segStartBase, segEndBase, wMedian / kCount);
                }
                segStartBase = contigs[i].start;
            }

            currentContig = contigs[i].contigName;
            segEndBase = contigs[i].end;
            wMedian = static_cast<float>(contigs[i].numRepeats) / coverage_depth;
            kCount = 1;
        }
    }

    if (segEndBase > segStartBase) {
        valSegment.emplace_back(currentContig, segStartBase, segEndBase, wMedian / kCount);
    }
    return valSegment;
}

void noiseFilter(CONTIGLIST& contigs, uint32_t windowSize)
{
    if (contigs.empty() || windowSize <= 1) {
        return;
    }

    const size_t cSize = contigs.size();
    const size_t half = windowSize / 2;
    std::vector<uint32_t> slidingWindow(windowSize);
    std::vector<uint32_t> filtered(cSize);

    for (size_t i = 0; i < cSize; i++)
    {
        const size_t winStart = (i >= half) ? (i - half) : 0;
        const size_t winEnd = std::min(cSize, i + (windowSize - half));
        const size_t wLen = winEnd - winStart;

        for (size_t j = 0; j < wLen; j++)
        {
            slidingWindow[j] = contigs[winStart + j].numRepeats;
        }

        std::sort(slidingWindow.begin(), slidingWindow.begin() + wLen);
        const uint32_t median = (wLen % 2 == 1)
            ? slidingWindow[wLen / 2]
            : (slidingWindow[wLen / 2 - 1] + slidingWindow[wLen / 2]) / 2;

        filtered[i] = median;
    }

    for (size_t i = 0; i < cSize; i++)
    {
        contigs[i].numRepeats = filtered[i];
    }
}

} // namespace

bool checkRepeats(
    const SEQLIST& assembly,
    const CBF* cbf,
    uint32_t kMerSize,
    const std::string& seed1,
    uint32_t windowSize,
    const std::string& outfile_path,
    uint32_t min_seg_len,
    uint32_t count_threshold,
    PipelineState& state)
{
    std::ofstream wout(outfile_path);
    if (!wout.is_open()) {
        std::cerr << "[!] Error: Could not open output BED file: " << outfile_path << std::endl;
        return false;
    }

    const float exact_coverage = (state.k_assembly > 0)
        ? (static_cast<float>(state.k_reads) / static_cast<float>(state.k_assembly))
        : 0.0f;
    const float safe_coverage = (state.k_assembly > 0 && state.k_reads > 0)
        ? exact_coverage
        : 1.0f;

    wout << "# k: " << kMerSize << " coverage depth: " << exact_coverage
         << " seg_min: " << min_seg_len << " gap_len: " << kMerSize
         << " count_threshold: " << count_threshold << "\n";

    const size_t num_elements = assembly.size();
    const std::vector<std::string> seeds = {seed1};
    const std::vector<btllib::SpacedSeed> parsedSeeds = btllib::parse_seeds(seeds);
    const unsigned cbfHashNum = cbf->get_hash_num();

    state.seg_table.clear();
    state.seg_table.resize(num_elements);

    #pragma omp parallel for schedule(dynamic, 1) shared(assembly, cbf, parsedSeeds, windowSize, kMerSize, count_threshold, safe_coverage, state)
    for (size_t idx = 0; idx < num_elements; idx++)
    {
        const std::string& contigSeq = assembly[idx].data;
        const std::string& contigName = assembly[idx].name;

        if (contigSeq.size() >= kMerSize)
        {
            CONTIGLIST contigs;
            btllib::SeedNtHash h(contigSeq, parsedSeeds, cbfHashNum, kMerSize);
            while (h.roll())
            {
                const auto count = cbf->contains(h.hashes());
                if (count != 0)
                {
                    const uint32_t start = static_cast<uint32_t>(h.get_pos());
                    const uint32_t end = start + kMerSize;
                    contigs.emplace_back(contigName, start, end, static_cast<uint32_t>(count));
                }
            }

            if (!contigs.empty())
            {
                noiseFilter(contigs, windowSize);
                state.seg_table[idx] = createSegments(contigs, kMerSize, count_threshold, safe_coverage);
            }
        }
    }

    for (size_t idx = 0; idx < num_elements; ++idx) {
        for (const SegmentInfo& seg : state.seg_table[idx]) {
            wout << seg.segmentName << "\t" << seg.start << "\t" << seg.end << "\t" << roundCount(seg.count) << "\n";
        }
    }

    wout.flush();
    const bool ok = wout.good();
    wout.close();
    return ok;
}
