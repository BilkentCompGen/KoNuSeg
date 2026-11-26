#include <checkRepeats.hpp>

// Regulate CPN to be inside user threshold
bool copyNumberInsideThreshold(float cn, float med, float th)
{	
	if (cn <= med + th && cn >= med - th) {
		return cn > 0.1f;
	}
	return false;
}

// Create unmerged segments for a single contig (lock-free; caller stores into pre-indexed segTable slot)
SEGMENTINFO* createSegments(
	CONTIGLIST *Contigs, uint k_len, uint count_threshold)	
{
	SEGMENTINFO *valSegment = new SEGMENTINFO();
	if (Contigs == nullptr || Contigs->empty()) {
		return valSegment;
	}

	string currentContig = Contigs->at(0).contigName;
	uint segStartBase = Contigs->at(0).start;
	uint segEndBase = Contigs->at(0).end;
	uint kCount = 1;
	const uint gap_len = k_len;
	const float coverage_depth = (k_assembly > 0 && k_reads > 0)
		? (static_cast<float>(k_reads) / static_cast<float>(k_assembly))
		: 1.0f;
	float wMedian = static_cast<float>(Contigs->at(0).numRepeats) / coverage_depth;

	// Iterate through all kmers in the contig (including the final kmer)
	for (size_t i = 1; i < Contigs->size(); i++)
	{	
		if (copyNumberInsideThreshold(static_cast<float>(Contigs->at(i).numRepeats) / coverage_depth, wMedian / kCount, count_threshold) && 
			(Contigs->at(i).contigName == currentContig) &&
			(Contigs->at(i).end >= segEndBase) &&
			(Contigs->at(i).end - segEndBase <= gap_len)
		) // Inside the threshold
		{		
			segEndBase = Contigs->at(i).end;
			wMedian += static_cast<float>(Contigs->at(i).numRepeats) / coverage_depth;
			kCount += 1;
		}
		else // Outside the threshold or gap too big 
		{
			if (Contigs->at(i).start < segEndBase && Contigs->at(i).contigName == currentContig)
			{
				if (segEndBase > segStartBase) {
					valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase, wMedian / kCount));
					segStartBase = segEndBase;
				}
			}
			else
			{
				if (segEndBase > segStartBase) {
					valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase, wMedian / kCount));
				}
				segStartBase = Contigs->at(i).start;
			}

			currentContig = Contigs->at(i).contigName;
			segEndBase = Contigs->at(i).end;
			wMedian = static_cast<float>(Contigs->at(i).numRepeats) / coverage_depth;
			kCount = 1;	
		}	
	} 	 

	// Always flush the final active segment of the contig using 0-based half-open [start, end) BED coordinates
	if (segEndBase > segStartBase) {
		valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase, wMedian / kCount));
	}
	return valSegment;
}


void noiseFilter(CONTIGLIST *contigs, uint windowSize)
{
	if (contigs == nullptr || contigs->empty() || windowSize <= 1) {
		return;
	}

	const size_t cSize = contigs->size();
	const size_t half = windowSize / 2;
	std::vector<uint> slidingWindow(windowSize);
	std::vector<uint> filtered(cSize);
	
	for (size_t i = 0; i < cSize; i++)
	{
		const size_t winStart = (i >= half) ? (i - half) : 0;
		const size_t winEnd = std::min(cSize, i + (windowSize - half));
		const size_t wLen = winEnd - winStart;

		for (size_t j = 0; j < wLen; j++)
		{
			slidingWindow[j] = contigs->at(winStart + j).numRepeats;
		}

		std::sort(slidingWindow.begin(), slidingWindow.begin() + wLen);
		uint median = (wLen % 2 == 1)
			? slidingWindow[wLen / 2]
			: (slidingWindow[wLen / 2 - 1] + slidingWindow[wLen / 2]) / 2;

		filtered[i] = median;
	}

	for (size_t i = 0; i < cSize; i++)
	{
		contigs->at(i).numRepeats = filtered[i];
	}
}


int query(const std::string& seq, const std::vector<std::string>& seeds, const btllib::CountingBloomFilter16 *cbf)
{
	if (seeds.empty() || cbf == nullptr) {
		return 0;
	}
	std::vector<int> repeats(seeds.size(), 0);
	btllib::SeedNtHash h(seq, seeds, cbf->get_hash_num(), seeds[0].size());

	while (h.roll()) {
		for (unsigned i = 0; i < seeds.size(); i++) {
			const auto count = cbf->contains(h.hashes() + i * cbf->get_hash_num());
			repeats[i] = static_cast<int>(count);
		}
	}
	return repeats[0];
}


bool checkRepeats(
	SEQLIST *assembly, CBF *cbf, uint kMerSize, string seed1, 
	uint windowSize, string outfile_path, uint min_seg_len, uint count_threshold)
{	
	ofstream wout(outfile_path);
	if (!wout.is_open()) {
		cerr << "[!] Error: Could not open output BED file: " << outfile_path << endl;
		return false;
	}

	const float exact_coverage = (k_assembly > 0)
		? (static_cast<float>(k_reads) / static_cast<float>(k_assembly))
		: 0.0f;
	wout << "# k: " << kMerSize << " coverage depth: " << exact_coverage
	     << " seg_min: " << min_seg_len << " gap_len: " << kMerSize
	     << " count_threshold: " << count_threshold << endl;

	const size_t num_elements = assembly->size();
	const std::vector<std::string> seeds = {seed1};
	const std::vector<btllib::SpacedSeed> parsedSeeds = btllib::parse_seeds(seeds);
	const unsigned cbfHashNum = cbf->get_hash_num();

	// Pre-size segTable so threads write deterministically without lock contention
	for (size_t i = 0; i < segTable->size(); ++i) {
		delete segTable->at(i);
	}
	segTable->clear();
	for (size_t i = 0; i < num_elements; ++i) {
		segTable->push_back(nullptr);
	}

	#pragma omp parallel for schedule(dynamic, 1) shared(assembly, cbf, parsedSeeds, windowSize, kMerSize, count_threshold)
	for (size_t idx = 0; idx < num_elements; idx++)
	{	
		const std::string& contigSeq = assembly->at(idx).data;
		const std::string& contigName = assembly->at(idx).name;

		if (contigSeq.size() >= kMerSize)
		{
			CONTIGLIST contigs;
			btllib::SeedNtHash h(contigSeq, parsedSeeds, cbfHashNum, kMerSize);
			while (h.roll())
			{
				const auto count = cbf->contains(h.hashes());
				if (count != 0)
				{
					const uint start = static_cast<uint>(h.get_pos());
					const uint end = start + kMerSize;
					contigs.push_back(ContigRepeatInfo(contigName, start, end, static_cast<uint>(count)));
				}
			}

			if (!contigs.empty())
			{
				noiseFilter(&contigs, windowSize);
				segTable->at(idx) = createSegments(&contigs, kMerSize, count_threshold);
			}
		}
	}

	// Write unmerged segments sequentially in deterministic FASTA contig order
	for (size_t idx = 0; idx < num_elements; ++idx) {
		SEGMENTINFO* segs = segTable->at(idx);
		if (segs == nullptr) continue;
		for (size_t s = 0; s < segs->size(); ++s) {
			const SegmentInfo& seg = segs->at(s);
			wout << seg.segmentName << "\t" << seg.start << "\t" << seg.end << "\t" << roundCount(seg.count) << "\n";
		}
	}

	wout.flush();
	bool ok = wout.good();
	wout.close();

	return ok;
}
