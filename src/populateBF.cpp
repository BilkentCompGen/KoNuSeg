#include <populateBF.hpp>

namespace {
struct SeqSlice {
	const char* data;
	size_t len;
};

constexpr size_t SLICE_KMER_TARGET = 65536;

std::vector<SeqSlice> buildAssemblySlices(SEQLIST *assembly, uint kMerSize) {
	std::vector<SeqSlice> slices;
	if (assembly == nullptr) return slices;
	slices.reserve(assembly->size());

	for (size_t idx = 0; idx < assembly->size(); idx++) {
		const std::string& seq = assembly->at(idx).data;
		if (seq.size() < kMerSize) continue;

		const size_t totalStarts = seq.size() - kMerSize + 1;
		for (size_t offset = 0; offset < totalStarts; offset += SLICE_KMER_TARGET) {
			const size_t startsInSlice = std::min(SLICE_KMER_TARGET, totalStarts - offset);
			const size_t sliceLen = startsInSlice + kMerSize - 1;
			slices.push_back(SeqSlice{seq.c_str() + offset, sliceLen});
		}
	}
	return slices;
}
} // namespace

ulong countAssemblyKmers(SEQLIST *assembly, uint kMerSize)
{
	std::vector<SeqSlice> slices = buildAssemblySlices(assembly, kMerSize);
	ulong totalValid = 0;

	#pragma omp parallel for schedule(dynamic, 1) reduction(+ : totalValid)
	for (size_t i = 0; i < slices.size(); i++) {
		btllib::NtHash h(slices[i].data, slices[i].len, 1, kMerSize);
		while (h.roll()) {
			totalValid += 1;
		}
	}
	return totalValid;
}

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

	// Create the SBF using bfBytes (btllib expects filter size in bytes), kMerSize, and optHashFunctions
	SBF *seed_bf = new SBF(bfBytes, kMerSize, {seed1}, optHashFunctions);

	// Populate Bloom filter using rolling SeedNtHash over fine-grained contig slices
	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >Loading Bloom filter..." << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	const std::vector<SeqSlice> slices = buildAssemblySlices(assembly, kMerSize);
	const auto& parsedSeeds = seed_bf->get_parsed_seeds();
	const unsigned hashNumPerSeed = seed_bf->get_hash_num_per_seed();
	const unsigned sbfK = seed_bf->get_k();
	const size_t numSeeds = parsedSeeds.size();

	#pragma omp parallel for schedule(dynamic, 1) shared(slices, seed_bf, parsedSeeds) reduction(+ : k_assembly)
	for (size_t i = 0; i < slices.size(); i++)
	{
		btllib::SeedNtHash h(slices[i].data, slices[i].len, parsedSeeds, hashNumPerSeed, sbfK);
		while (h.roll())
		{
			for (size_t s = 0; s < numSeeds; s++)
			{
				seed_bf->insert(h.hashes() + s * hashNumPerSeed);
			}
			k_assembly += 1;
		}
	}

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >BF occupancy: " << seed_bf->get_occupancy() << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;
	return seed_bf;
}
