#include <loadReads.hpp>


CBF *loadReads(
	string readsFile, uint k_len, SEQLIST *assembly, SBF *seed_bf, 	ulong k_num, double fpr, string seed1,uint kMerSize)
{
	// Variables for reading sequences
	string line, name, content; 
	uint read_counter = 0;

	
	// Calculate the optimal size (m)
	size_t cbfSize = static_cast<size_t>((-1) * (k_num * log(fpr)) / (log(2) * log(2)));

	// Calculate and store the optimal number of hash functions (k)
	uint64_t optHashFunctions = static_cast<uint64_t>(round(((double)cbfSize / k_num) * log(2)));

	// Safety check: ensure at least one hash function is used
	if (optHashFunctions < 1) optHashFunctions = 1;

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >CBF Size:    " << cbfSize <<  " bits" << endl;

	// Calculate GB for printing: CBF uses 1 byte per slot
	cout << "\t\t\t[i] >CBF Memory:   " << std::fixed << std::setprecision(4) 
		<< (double)cbfSize / (1024.0 * 1024.0 * 1024.0) << " GB" << endl;

	cout << "\t\t\t[i] >Optimal Hash Functions: " << optHashFunctions << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	// Initialize the CBF with the optimized values
	CBF* cbf = new CBF(cbfSize, optHashFunctions);

	ifstream reads(readsFile);

	
	// Get Reads
	while (std::getline(reads, line)) {
		if (line.empty()) continue;

		if (line[0] == '@') {
			// If we have DNA from the previous read, save it now
			if (!content.empty()) {
				populateCBF(cbf, content, assembly, seed_bf, k_num, fpr, seed1, kMerSize);
				content = ""; // Clear for the next read
			}
			
			// Start the new read
			boost::trim_right(line);
			name = line;
			read_counter = 1; 
		}
		else if (line[0] == '+') {
			read_counter = 0; // Stop collecting DNA
		}
		else if (read_counter == 1) {
			// STEP 3: Accumulate multi-line DNA
			boost::trim_right(line);
			content.append(line);
		}
	}

	// Process the very last read
	if (!content.empty()) {
		populateCBF(cbf, content, assembly, seed_bf, k_num, fpr, seed1, kMerSize);
	}

	cout <<"\t\t\t[i] >----------------------------------------" << endl;
	cout <<"\t\t\t[i] >CBF occupancy: " << cbf->get_occupancy() << endl;
	cout <<"\t\t\t[i] >----------------------------------------" << endl;

	cout <<"\t\t\t[i] >----------------------------------------" << endl;
	cout <<"\t\t\t[i] >K-MERS IN BLOOM FILTER: " << k_assembly << endl;
	cout <<"\t\t\t[i] >K-MERS IN COUNTING BLOOM FILTER: " << k_reads << endl;
	cout <<"\t\t\t[i] >BLOOM FITLER COVERAGE DEPTH: " << round(static_cast<float>(k_reads) / k_assembly) << endl;
	cout <<"\t\t\t[i] >----------------------------------------" << endl;

	return cbf;
}




void populateCBF(
	CBF *cbf, string read, SEQLIST *assembly, SBF *seed_bf, ulong k_num, double fpr, string seed1,uint kMerSize)
{
	string temp = "";
  	std::vector<std::string> seeds = {seed1};
	if (read.size() >= kMerSize)
	{
		#pragma omp parallel for shared(read, cbf, seed_bf) private(temp) reduction(+:k_reads)
		for (uint j = 0; j <= read.size() - kMerSize; j++)
		{
			temp = read.substr(j, kMerSize);

			auto hit_seeds = seed_bf->contains(temp);

			if (find(hit_seeds[0].begin(), hit_seeds[0].end(), 0) != hit_seeds[0].end())
			{
				populate(temp, seeds, cbf);
				k_reads = k_reads + 1;
			}
		}
	}
}


void populate(const std::string& seq, const std::vector<std::string>& seeds, btllib::CountingBloomFilter16 *cbf)
{
	btllib::SeedNtHash h(seq, seeds, cbf->get_hash_num(), seeds[0].size());
	while (h.roll()) {
		for (unsigned i = 0; i < seeds.size(); i++) {
		// temp array for seed[i]'s hashes
			std::unique_ptr<uint64_t[]> hashes(new uint64_t[cbf->get_hash_num()]);
			std::copy(h.hashes() + i * cbf->get_hash_num(), h.hashes() + (i + 1) * cbf->get_hash_num(), hashes.get());
			cbf->insert(hashes.get());
		}
	}
}

