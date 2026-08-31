#include <populateBF.hpp>



SBF *populateBF(SEQLIST *assembly, string seed1, ulong a_num, double fpr, uint kMerSize)
{

	uint s_pos = 0;

	// Calculate the optimal size (m)
	size_t bfSize = static_cast<size_t>((-1) * (a_num * log(fpr)) / (log(2) * log(2)));

	// Calculate the optimal number of hash functions (k). We use round() because k must be an integer
	uint64_t optHashFunctions = static_cast<uint64_t>(round((static_cast<double>(bfSize) / a_num) * log(2)));

	// Ensure at least 1 hash function is used
	if (optHashFunctions < 1) optHashFunctions = 1;

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >BF size: " << bfSize << " bits" << endl;
	cout << "\t\t\t[i] >Optimal Hash Functions: " << optHashFunctions << endl;
	cout << "\t\t\t[i] >Standard BF Memory: " << (double)bfSize / (8.0 * 1024 * 1024 * 1024) << " GB" << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	// Create the SBF using the calculated optHashFunctions
	SBF *seed_bf = new SBF(bfSize, seed1.size(), {seed1}, optHashFunctions);

	// Populate Bloom filter
	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >Loading Bloom filter..." << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	uint num_elements = assembly->size();
	uint idx;

	for (idx = 0; idx < num_elements; idx++)
	{
		if (assembly->at(idx).data.size() >= kMerSize)
		{
			uint j;
     		#pragma omp parallel for shared(assembly, seed_bf) private(j) reduction(+ : k_assembly)
			for (j = 0; j <= assembly->at(idx).data.size() - kMerSize; j++)
			{
				seed_bf->insert(assembly->at(idx).data.substr(j, kMerSize));				
				s_pos = j + kMerSize;
				k_assembly = k_assembly + 1;
			}
		}

	}

	cout <<"\t\t\t[i] >----------------------------------------" << endl;
	cout <<"\t\t\t[i] >BF occupancy: " << seed_bf->get_occupancy() << endl;
	cout <<"\t\t\t[i] >----------------------------------------" << endl;
	return seed_bf;
}
