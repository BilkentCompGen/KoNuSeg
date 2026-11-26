#include <populateBF.hpp>

namespace {
struct SeqSlice {
    const char* data;
    size_t len;
};

constexpr size_t SLICE_KMER_TARGET = 65536;

std::vector<SeqSlice> buildAssemblySlices(const SEQLIST& assembly, uint32_t kMerSize) {
    std::vector<SeqSlice> slices;
    slices.reserve(assembly.size());

    for (const Sequence& rec : assembly) {
        const std::string& seq = rec.data;
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

uint64_t countAssemblyKmers(const SEQLIST& assembly, uint32_t kMerSize)
{
    const std::vector<SeqSlice> slices = buildAssemblySlices(assembly, kMerSize);
    uint64_t totalValid = 0;

    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : totalValid)
    for (size_t i = 0; i < slices.size(); i++) {
        btllib::NtHash h(slices[i].data, slices[i].len, 1, kMerSize);
        while (h.roll()) {
            totalValid += 1;
        }
    }
    return totalValid;
}

std::unique_ptr<SBF> populateBF(
    const SEQLIST& assembly,
    const std::string& seed1,
    uint64_t a_num,
    double fpr,
    uint32_t kMerSize,
    uint64_t& k_assembly_out)
{
    const double ln2 = std::log(2.0);
    size_t bfBits = static_cast<size_t>((-1.0) * (static_cast<double>(a_num) * std::log(fpr)) / (ln2 * ln2));
    if (bfBits < CHAR_BIT) bfBits = CHAR_BIT;
    const size_t bfBytes = (bfBits + CHAR_BIT - 1) / CHAR_BIT;

    uint64_t optHashFunctions = static_cast<uint64_t>(std::round((static_cast<double>(bfBits) / static_cast<double>(a_num)) * ln2));
    if (optHashFunctions < 1) optHashFunctions = 1;

    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    std::cout << "\t\t\t[i] >BF size: " << bfBits << " bits (" << bfBytes << " bytes)" << std::endl;
    std::cout << "\t\t\t[i] >Optimal Hash Functions: " << optHashFunctions << std::endl;
    std::cout << "\t\t\t[i] >Standard BF Memory: " << static_cast<double>(bfBytes) / (1024.0 * 1024.0 * 1024.0) << " GB" << std::endl;
    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;

    auto seed_bf = std::make_unique<SBF>(bfBytes, kMerSize, std::vector<std::string>{seed1}, static_cast<unsigned>(optHashFunctions));

    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    std::cout << "\t\t\t[i] >Loading Bloom filter..." << std::endl;
    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;

    const std::vector<SeqSlice> slices = buildAssemblySlices(assembly, kMerSize);
    const auto& parsedSeeds = seed_bf->get_parsed_seeds();
    const unsigned hashNumPerSeed = seed_bf->get_hash_num_per_seed();
    const unsigned sbfK = seed_bf->get_k();
    const size_t numSeeds = parsedSeeds.size();
    uint64_t totalKmers = 0;

    #pragma omp parallel for schedule(dynamic, 1) shared(slices, seed_bf, parsedSeeds) reduction(+ : totalKmers)
    for (size_t i = 0; i < slices.size(); i++)
    {
        btllib::SeedNtHash h(slices[i].data, slices[i].len, parsedSeeds, hashNumPerSeed, sbfK);
        while (h.roll())
        {
            for (size_t s = 0; s < numSeeds; s++)
            {
                seed_bf->insert(h.hashes() + s * hashNumPerSeed);
            }
            totalKmers += 1;
        }
    }

    k_assembly_out = totalKmers;

    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    std::cout << "\t\t\t[i] >BF occupancy: " << seed_bf->get_occupancy() << std::endl;
    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    return seed_bf;
}
