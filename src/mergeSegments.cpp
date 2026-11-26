#include <mergeSegments.hpp>

namespace {

// Sort pooled copy estimates, discard the top quartile (at least 1 highest outlier when size >= 2)
// to suppress Bloom filter hash-collision inflation, and return the rounded median of the remainder.
float processData(std::vector<float>& arr) {
    if (arr.empty()) {
        std::cerr << "Error: No values to calculate median" << std::endl;
        return 0.0f;
    }

    std::sort(arr.begin(), arr.end());

    const size_t size = arr.size();
    const size_t remove_count = (size >= 2) ? ((size + 3) / 4) : 0;
    const size_t remaining_size = size - remove_count;

    if (remaining_size % 2 == 0) {
        return roundCount((arr[remaining_size / 2 - 1] + arr[remaining_size / 2]) / 2.0f);
    }
    return roundCount(arr[remaining_size / 2]);
}

} // namespace

bool mergeSegments(
    const std::string& outfile_path,
    uint32_t k_len,
    uint32_t min_seg_len,
    uint32_t count_threshold,
    const SEGMENTTABLE& seg_table)
{
    const std::string merged_path = konuseg::deriveMergedBedPath(outfile_path);
    std::ofstream wout(merged_path);
    if (!wout.is_open()) {
        std::cerr << "[!] Error: Could not open merged output BED file: " << merged_path << std::endl;
        return false;
    }

    std::string currentSegment;
    uint32_t segStartBase = 0;
    uint32_t segEndBase = 0;
    std::vector<float> values;
    SEGMENTINFO valSegment;

    for (const SEGMENTINFO& segItem : seg_table)
    {
        for (const SegmentInfo& nextSeg : segItem)
        {
            const float nextRounded = roundCount(nextSeg.count);

            if (currentSegment.empty())
            {
                currentSegment = nextSeg.segmentName;
                segStartBase = nextSeg.start;
                segEndBase = nextSeg.end;
                values.push_back(nextRounded);
            }
            else if (currentSegment == nextSeg.segmentName)
            {
                const uint32_t currentLen = (segEndBase >= segStartBase) ? (segEndBase - segStartBase) : 0;
                const uint32_t nextLen = (nextSeg.end >= nextSeg.start) ? (nextSeg.end - nextSeg.start) : 0;
                const bool withinGap = (nextSeg.start <= segEndBase) ||
                    ((nextSeg.start - segEndBase) <= 3 * k_len);

                // Avoid absorbing a >= k_len flanking locus into an already full-sized (>= min_seg_len) block
                // when their copy numbers differ by more than count_threshold, while smoothing sub-kmer (< k_len) artifacts.
                const bool incompatibleLongBlock = (currentLen >= k_len) &&
                    (nextLen >= min_seg_len) &&
                    !values.empty() &&
                    (std::fabs(nextRounded - values.back()) > static_cast<float>(count_threshold));

                if (currentLen < min_seg_len && withinGap && !incompatibleLongBlock)
                {
                    if (nextSeg.end > segEndBase) {
                        segEndBase = nextSeg.end;
                    }
                    values.push_back(nextRounded);
                }
                else
                {
                    const float mergedCount = processData(values);
                    if (mergedCount > 0.0f) {
                        valSegment.emplace_back(currentSegment, segStartBase, segEndBase, mergedCount);
                    }

                    currentSegment = nextSeg.segmentName;
                    segStartBase = nextSeg.start;
                    segEndBase = nextSeg.end;
                    values.clear();
                    values.push_back(nextRounded);
                }
            }
            else
            {
                const float mergedCount = processData(values);
                if (mergedCount > 0.0f) {
                    valSegment.emplace_back(currentSegment, segStartBase, segEndBase, mergedCount);
                }
                currentSegment = nextSeg.segmentName;
                segStartBase = nextSeg.start;
                segEndBase = nextSeg.end;
                values.clear();
                values.push_back(nextRounded);
            }
        }
    }

    if (!currentSegment.empty() && !values.empty())
    {
        const float mergedCount = processData(values);
        if (mergedCount > 0.0f) {
            valSegment.emplace_back(currentSegment, segStartBase, segEndBase, mergedCount);
        }
    }

    for (const SegmentInfo& seg : valSegment)
    {
        wout << seg.segmentName << "\t" << seg.start << "\t" << seg.end << "\t" << seg.count << "\n";
    }

    wout.flush();
    const bool ok = wout.good();
    wout.close();
    return ok;
}
