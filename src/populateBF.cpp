#include <populateBF.hpp>



SBF *populateBF(SEQLIST *assembly, string seed1, ulong a_num, double fpr, uint kMerSize)
{
	// Calculate the optimal size in bits (m) and convert to bytes for btllib::SeedBloomFilter
	size_t bfBits = static_cast<size_t>((-1.0) * (static_cast<double>(a_num) * log(fpr)) / (log(2) * log(2)));
	if (bfBits < CHAR_BIT) bfBits = CHAR_BIT;
	size_t bfBytes = (bfBits + CHAR_BIT - 1) / CHAR_BIT;

	// Calculate the optimal number of hash functions (k). We use round() because k must be an integer
	uint64_t optHashFunctions = static_cast<uint64_t>(round((static_cast<double>(bfBits) / a_num) * log(2)));

	// Ensure at least 1 hash function is used
	if (optHashFunctions < 1) optHashFunctions = 1;

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >BF size: " << bfBits << " bits (" << bfBytes << " bytes)" << endl;
	cout << "\t\t\t[i] >Optimal Hash Functions: " << optHashFunctions << endl;
	cout << "\t\t\t[i] >Standard BF Memory: " << static_cast<double>(bfBytes) / (1024.0 * 1024.0 * 1024.0) << " GB" << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	// Create the SBF using bfBytes (btllib expects filter size in bytes) and optHashFunctions
	SBF *seed_bf = new SBF(bfBytes, seed1.size(), {seed1}, optHashFunctions);

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
				string kmer = assembly->at(idx).data.substr(j, kMerSize);
				btllib::SeedNtHash h(kmer, seed_bf->get_parsed_seeds(), seed_bf->get_hash_num_per_seed(), seed_bf->get_k());
				if (h.roll())
				{
					for (size_t s = 0; s < seed_bf->get_parsed_seeds().size(); s++)
					{
						seed_bf->insert(h.hashes() + s * seed_bf->get_hash_num_per_seed());
					}
					k_assembly = k_assembly + 1;
				}
			}
		}

	}

	cout <<"\t\t\t[i] >----------------------------------------" << endl;
	cout <<"\t\t\t[i] >BF occupancy: " << seed_bf->get_occupancy() << endl;
	cout <<"\t\t\t[i] >----------------------------------------" << endl;
	return seed_bf;
}
