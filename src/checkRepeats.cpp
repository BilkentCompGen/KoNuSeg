#include <checkRepeats.hpp>

// Regulate CPN to be inside user threshold
bool copyNumberInsideThreshold(float cn, float med, float th)
{	
	if(cn <= med + th && cn >= med - th ){
		if(cn > 0.1){
			return true;
		}
		else{
			return false;		
		}
	}
	else{
		return false;
	}

}

// Create unmerged segments
void createSegments(
	CONTIGLIST *Contigs, ofstream* wout, uint k_len, uint min_seg_len, uint count_threshold)	
{
	if (Contigs == nullptr || Contigs->empty()) {
		return;
	}

	#pragma omp critical
	{	
		string currentContig = "";
		uint segStartBase, segEndBase;
		float wMedian;
		uint kCount;
		uint gap_len = k_len;
		const float coverage_depth = (k_assembly > 0 && k_reads > 0)
			? (static_cast<float>(k_reads) / static_cast<float>(k_assembly))
			: 1.0f;
		
		// Initialising variables
		currentContig = Contigs->at(0).contigName;
		segStartBase = Contigs->at(0).start;
		segEndBase = Contigs->at(0).end;
		kCount = 1;
		wMedian = static_cast<float>(Contigs->at(0).numRepeats) / coverage_depth;

		SEGMENTINFO *valSegment = new SEGMENTINFO();

		// Iterating through all kmers in the contig (including the final kmer)
		for (size_t i = 1; i < Contigs->size(); i++)
		{	
			if (copyNumberInsideThreshold(static_cast<float>(Contigs->at(i).numRepeats) / coverage_depth, wMedian / kCount, count_threshold) && 
				(Contigs->at(i).contigName == currentContig) &&
				(Contigs->at(i).end >= segEndBase) &&
				(Contigs->at(i).end - segEndBase <= gap_len)
			) // Inside the threshold
			{		
				segEndBase = Contigs->at(i).end;
				
				wMedian = wMedian + (static_cast<float>(Contigs->at(i).numRepeats) / coverage_depth);
				kCount = kCount + 1;
			}
			else // Outside the threshold or gap too big 
			{
				if (Contigs->at(i).start < segEndBase && Contigs->at(i).contigName == currentContig)
				{
					if (segEndBase > segStartBase + 1) {
						valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase - 1, wMedian / kCount));
						*wout << currentContig << "\t" << segStartBase << "\t" << segEndBase - 1 << "\t" << roundCount(wMedian / kCount) << endl;
						segStartBase = segEndBase;
					}
				}
				else
				{
					if (segEndBase > segStartBase + 1) {
						valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase - 1, wMedian / kCount));
						*wout << currentContig << "\t" << segStartBase << "\t" << segEndBase - 1 << "\t" << roundCount(wMedian / kCount) << endl;
					}
					segStartBase = Contigs->at(i).start;
				}

				currentContig = Contigs->at(i).contigName;
				segEndBase = Contigs->at(i).end;
				
				wMedian = static_cast<float>(Contigs->at(i).numRepeats) / coverage_depth;
				kCount = 1;	
			}	
		} 	 

		// Always flush the final active segment of the contig
		if (segEndBase > segStartBase + 1) {
			valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase - 1, wMedian / kCount));
			*wout << currentContig << "\t" << segStartBase << "\t" << segEndBase - 1 << "\t" << roundCount(wMedian / kCount) << endl;
		}
		segTable->push_back(valSegment);
	}
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
	if (seeds.empty()) {
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






void checkRepeats(
	SEQLIST *assembly, CBF *cbf, uint kMerSize, string seed1, 
	long double genomeSize, uint windowSize, string outfile_path, uint min_seg_len, uint count_threshold)
{	
	// Variables for reading sequences
	string name = "";
	string temp = "";

	uint repeats = 0;
	uint start = 0;
	uint end = 0;
	uint num_elements = assembly->size();
	uint idx;
	std::vector<std::string> seeds = {seed1};

	ofstream wout(outfile_path);
	const float header_coverage = (k_assembly > 0)
		? (static_cast<float>(k_reads) / static_cast<float>(k_assembly))
		: 0.0f;
	wout << "# k: " << kMerSize << " coverage depth: " << roundCount(header_coverage) << " seg_min: " << min_seg_len << " gap_len: " << kMerSize << " count_threshold: " << count_threshold << endl;


	#pragma omp parallel for shared(assembly, cbf, windowSize, outfile_path, kMerSize, min_seg_len, count_threshold) private (idx, repeats, temp, start, end)
	for (idx = 0; idx < num_elements; idx++)
	{	
		CONTIGLIST *contigs = new CONTIGLIST();
		if (assembly->at(idx).data.size() >= kMerSize)
		{
			for (uint j = 0; j <= assembly->at(idx).data.size() - kMerSize; j++)
			{
				temp = assembly->at(idx).data.substr(j, kMerSize);
				repeats = query(temp, seeds, cbf);
				if (repeats != 0 )
				{
					start = j;
					end = start + kMerSize;
					contigs->push_back(ContigRepeatInfo(assembly->at(idx).name.substr(1, assembly->at(idx).name.size()), start, end, repeats));
				}
			}

			if (contigs->size() > 0)
			{
				// First filter out the noise in the contig
				noiseFilter(contigs, windowSize);
				// Create Segments in the contig
				createSegments(contigs, &wout, kMerSize, min_seg_len, count_threshold);
			}
			
		} 
		delete contigs;

 	} 	

	 wout.close();

	return;

}
