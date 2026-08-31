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
	#pragma omp critical
	{	
	
		string currentContig = "";
		uint segStartBase, segEndBase, currentBase;
		float wMedian;
		uint kCount;
		uint gap_len = k_len;
		uint first_check=1;
		uint coverage_depth = roundCount(static_cast<float>(k_reads) / k_assembly);
		
		// Initialising variables
		currentContig = Contigs->at(0).contigName;
		segStartBase = Contigs->at(0).start;
		segEndBase = Contigs->at(0).end;
		kCount = 1;
		wMedian = Contigs->at(0).numRepeats/coverage_depth;

		SEGMENTINFO *valSegment= new SEGMENTINFO();

		// Iterating through kmers
		for (uint i = 1; i <= Contigs->size(); i++)
		{	
			if (i< Contigs->size()-1)
			{	

				if(copyNumberInsideThreshold(Contigs->at(i).numRepeats / coverage_depth, wMedian/kCount, count_threshold) && 
					(Contigs->at(i).end-segEndBase <= gap_len)
				) // Inside the threshold
				{		
						currentBase = Contigs->at(i).start;
						segEndBase = Contigs->at(i).end;
						
						wMedian = wMedian + (Contigs->at(i).numRepeats * 1.0 / coverage_depth * 1.0);
						kCount = kCount + 1;
				}
				else // Outside the threshold or gap too big 
				{
					first_check=0;
					if (Contigs->at(i).start < segEndBase && Contigs->at(i).contigName.compare(Contigs->at(i).contigName) == 0 )
					{
						if( (Contigs->at(i+1).start - 2) > segStartBase){
							valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase-1 , wMedian/kCount));
							*wout<< currentContig << "\t" << segStartBase << "\t" << segEndBase -1 << "\t" << roundCount(wMedian/kCount) << endl;
							segStartBase = segEndBase;
						}											
					}
					else
					{
						if(segEndBase>segStartBase){
							valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase-1, wMedian/kCount));
							*wout<< currentContig << "\t" << segStartBase << "\t" << segEndBase -1 << "\t" << roundCount(wMedian/kCount) << endl;
							segStartBase = Contigs->at(i+1).start;
						}
					}

					currentContig = Contigs->at(i).contigName;
					currentBase = Contigs->at(i).start;
					segEndBase = Contigs->at(i).end;
					
					wMedian = (Contigs->at(i).numRepeats * 1.0 / coverage_depth * 1.0);
					kCount = 1;	
				}	
	
			}
		} 	 

		if(first_check==1){
			valSegment->push_back(SegmentInfo(currentContig, segStartBase, segEndBase-1, wMedian/kCount));
			*wout<< currentContig << "\t" << segStartBase << "\t" << segEndBase -1 << "\t" << roundCount(wMedian/kCount) << endl;
		}
		segTable->push_back(valSegment);
	}
}



void noiseFilter(CONTIGLIST *contigs, uint windowSize)
{
	uint wIndex;
	uint cSize = contigs->size();
	uint slidingWindow[windowSize];
	
	for (uint i = 0; i < cSize; i++)
	{
		wIndex = 0;
		for (uint j = 0; j <= windowSize; j++)
		{
			if (i + j < cSize)
			{
				slidingWindow[j] = contigs->at(i + j).numRepeats;
				wIndex++;
			}
		}

		/*Here we take two parameters, the beginning of the
		array and the length n up to which we want the array to
		be sorted*/
		if (wIndex == windowSize )
		{
			sort(slidingWindow, slidingWindow + windowSize);
			// Find median of sliding window
			uint median = slidingWindow[(windowSize / 2)];
	
			// Assign median value of window to mid index
			contigs->at(i + (windowSize / 2)).numRepeats = median;
		}
	}
	return;
}




int query(const std::string& seq, const std::vector<std::string>& seeds, const btllib::CountingBloomFilter16 *cbf)
{
	int repeats[seeds.size()];
	btllib::SeedNtHash h(seq, seeds, cbf->get_hash_num(), seeds[0].size());

		while (h.roll()) {
			for (unsigned i = 0; i < seeds.size(); i++) {
				std::unique_ptr<uint64_t[]> hashes(new uint64_t[cbf->get_hash_num()]);
				std::copy(h.hashes() + i * cbf->get_hash_num(), h.hashes() + (i + 1) * cbf->get_hash_num(), hashes.get());
				const auto count = cbf->contains(hashes.get());
				repeats[i]=(unsigned)count;
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
	wout << "# k: " << kMerSize<< " coverage depth: " << roundCount(static_cast<float>(k_reads) / k_assembly) << " seg_min: " << min_seg_len << " gap_len: " << kMerSize<< " count_threshold: " << count_threshold << endl;


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
