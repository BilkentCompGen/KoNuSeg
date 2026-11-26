#include <loadReads.hpp>

namespace {
struct ReadSlice {
    const char* data;
    size_t len;
};

constexpr size_t READ_SLICE_KMER_TARGET = 65536;
constexpr size_t BATCH_MAX_BYTES = 64ULL * 1024ULL * 1024ULL; // 64 MiB max buffered reads per batch

uint64_t processReadSlice(
    const char* data,
    size_t len,
    CBF* cbf,
    const SBF* seed_bf,
    const std::vector<btllib::SpacedSeed>& parsedCbfSeeds,
    unsigned cbfHashNum,
    unsigned sbfHashNum,
    bool reuseSbfHashes,
    uint32_t kMerSize)
{
    uint64_t hits = 0;
    if (len < kMerSize) return 0;

    if (reuseSbfHashes) {
        btllib::SeedNtHash h(data, len, seed_bf->get_parsed_seeds(), sbfHashNum, kMerSize);
        while (h.roll()) {
            if (seed_bf->contains(h.hashes())) {
                cbf->insert(h.hashes());
                hits += 1;
            }
        }
    } else {
        btllib::SeedNtHash hSbf(data, len, seed_bf->get_parsed_seeds(), sbfHashNum, seed_bf->get_k());
        while (hSbf.roll()) {
            if (seed_bf->contains(hSbf.hashes())) {
                const size_t kmerPos = hSbf.get_pos();
                btllib::SeedNtHash hCbf(data + kmerPos, kMerSize, parsedCbfSeeds, cbfHashNum, kMerSize);
                if (hCbf.roll()) {
                    cbf->insert(hCbf.hashes());
                    hits += 1;
                }
            }
        }
    }
    return hits;
}

uint64_t flushReadBatch(
    const std::vector<std::string>& batch,
    CBF* cbf,
    const SBF* seed_bf,
    const std::vector<btllib::SpacedSeed>& parsedCbfSeeds,
    bool reuseSbfHashes,
    uint32_t kMerSize)
{
    if (batch.empty()) return 0;

    std::vector<ReadSlice> slices;
    slices.reserve(batch.size());
    for (const std::string& r : batch) {
        if (r.size() < kMerSize) continue;
        const size_t totalStarts = r.size() - kMerSize + 1;
        for (size_t offset = 0; offset < totalStarts; offset += READ_SLICE_KMER_TARGET) {
            const size_t startsInSlice = std::min(READ_SLICE_KMER_TARGET, totalStarts - offset);
            const size_t sliceLen = startsInSlice + kMerSize - 1;
            slices.push_back(ReadSlice{r.c_str() + offset, sliceLen});
        }
    }

    const unsigned cbfHashNum = cbf->get_hash_num();
    const unsigned sbfHashNum = seed_bf->get_hash_num_per_seed();
    uint64_t batchHits = 0;

    #pragma omp parallel for schedule(dynamic, 1) shared(slices, cbf, seed_bf, parsedCbfSeeds) reduction(+ : batchHits)
    for (size_t i = 0; i < slices.size(); i++) {
        batchHits += processReadSlice(
            slices[i].data,
            slices[i].len,
            cbf,
            seed_bf,
            parsedCbfSeeds,
            cbfHashNum,
            sbfHashNum,
            reuseSbfHashes,
            kMerSize);
    }
    return batchHits;
}
} // namespace

std::unique_ptr<CBF> loadReads(
    const std::string& readsFile,
    const SBF* seed_bf,
    uint64_t k_num,
    double fpr,
    const std::string& seed1,
    uint32_t kMerSize,
    uint64_t k_assembly,
    uint64_t& k_reads_out)
{
    if (readsFile != "-" && readsFile.find("://") == std::string::npos) {
        std::ifstream checkFile(readsFile);
        if (!checkFile.is_open()) {
            std::cerr << "[!] Error: Could not open reads file: " << readsFile << std::endl;
            return nullptr;
        }
    }

    btllib::DataSource readsSource(readsFile);
    FILE* fp = readsSource;
    if (fp == nullptr) {
        std::cerr << "[!] Error: Could not open reads stream: " << readsFile << std::endl;
        return nullptr;
    }

    const double ln2 = std::log(2.0);
    size_t cbfSlots = static_cast<size_t>((-1.0) * (static_cast<double>(k_num) * std::log(fpr)) / (ln2 * ln2));
    if (cbfSlots < 1) cbfSlots = 1;
    const size_t cbfBytes = cbfSlots * sizeof(uint16_t);

    uint64_t optHashFunctions = static_cast<uint64_t>(std::round((static_cast<double>(cbfSlots) / static_cast<double>(k_num)) * ln2));
    if (optHashFunctions < 1) optHashFunctions = 1;

    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    std::cout << "\t\t\t[i] >CBF Size:    " << cbfSlots << " counters (" << cbfBytes << " bytes)" << std::endl;
    std::cout << "\t\t\t[i] >CBF Memory:   " << std::fixed << std::setprecision(4)
              << static_cast<double>(cbfBytes) / (1024.0 * 1024.0 * 1024.0) << " GB" << std::endl;
    std::cout << "\t\t\t[i] >Optimal Hash Functions: " << optHashFunctions << std::endl;
    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;

    auto cbf = std::make_unique<CBF>(cbfBytes, static_cast<unsigned>(optHashFunctions));

    const std::vector<std::string> seeds = {seed1};
    const std::vector<btllib::SpacedSeed> parsedCbfSeeds = btllib::parse_seeds(seeds);
    const bool reuseSbfHashes = (seed_bf != nullptr) &&
        (seed_bf->get_hash_num_per_seed() == cbf->get_hash_num()) &&
        (seed_bf->get_k() == kMerSize) &&
        (seed_bf->get_seeds() == seeds);

    const size_t batchMaxReads = std::max<size_t>(64, static_cast<size_t>(omp_get_max_threads()) * 8);
    std::vector<std::string> readBatch;
    readBatch.reserve(batchMaxReads);
    size_t batchBytes = 0;
    uint64_t totalReadKmers = 0;

    auto enqueueRead = [&](std::string& seq) {
        if (seq.empty()) return;
        batchBytes += seq.size();
        readBatch.push_back(std::move(seq));
        seq.clear();
        if (readBatch.size() >= batchMaxReads || batchBytes >= BATCH_MAX_BYTES) {
            totalReadKmers += flushReadBatch(readBatch, cbf.get(), seed_bf, parsedCbfSeeds, reuseSbfHashes, kMerSize);
            readBatch.clear();
            batchBytes = 0;
        }
    };

    enum class ReadState { HEADER, FASTQ_SEQ, FASTQ_QUAL, FASTA_SEQ };
    ReadState state = ReadState::HEADER;
    std::string line, content;
    size_t seqLines = 0;
    size_t qualLen = 0;

    while (konuseg::readLineFromStream(fp, line)) {
        konuseg::trimRight(line);
        if (line.empty()) continue;

        if (state == ReadState::HEADER) {
            if (line[0] == '@') {
                content.clear();
                seqLines = 0;
                qualLen = 0;
                state = ReadState::FASTQ_SEQ;
            } else if (line[0] == '>') {
                content.clear();
                state = ReadState::FASTA_SEQ;
            }
        } else if (state == ReadState::FASTQ_SEQ) {
            if (line[0] == '+') {
                qualLen = 0;
                state = ReadState::FASTQ_QUAL;
            } else {
                content.append(line);
                seqLines++;
            }
        } else if (state == ReadState::FASTQ_QUAL) {
            qualLen += line.size();
            if (qualLen >= content.size() || seqLines <= 1) {
                enqueueRead(content);
                state = ReadState::HEADER;
            }
        } else if (state == ReadState::FASTA_SEQ) {
            if (line[0] == '>') {
                enqueueRead(content);
            } else {
                content.append(line);
            }
        }
    }

    if (!content.empty() && (state == ReadState::FASTA_SEQ || state == ReadState::FASTQ_SEQ || state == ReadState::FASTQ_QUAL)) {
        enqueueRead(content);
    }
    if (!readBatch.empty()) {
        totalReadKmers += flushReadBatch(readBatch, cbf.get(), seed_bf, parsedCbfSeeds, reuseSbfHashes, kMerSize);
        readBatch.clear();
    }

    readsSource.close();
    k_reads_out = totalReadKmers;

    const double exactCoverage = (k_assembly > 0)
        ? (static_cast<double>(totalReadKmers) / static_cast<double>(k_assembly))
        : 0.0;

    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    std::cout << "\t\t\t[i] >CBF occupancy: " << cbf->get_occupancy() << std::endl;
    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;

    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;
    std::cout << "\t\t\t[i] >K-MERS IN BLOOM FILTER: " << k_assembly << std::endl;
    std::cout << "\t\t\t[i] >K-MERS IN COUNTING BLOOM FILTER: " << totalReadKmers << std::endl;
    std::cout << "\t\t\t[i] >BLOOM FILTER COVERAGE DEPTH: " << exactCoverage << std::endl;
    std::cout << "\t\t\t[i] >----------------------------------------" << std::endl;

    return cbf;
}

uint64_t populateCBF(
    CBF* cbf,
    const std::string& read,
    const SBF* seed_bf,
    const std::string& seed1,
    uint32_t kMerSize)
{
    if (read.size() < kMerSize || cbf == nullptr || seed_bf == nullptr) return 0;
    const std::vector<std::string> seeds = {seed1};
    const std::vector<btllib::SpacedSeed> parsedCbfSeeds = btllib::parse_seeds(seeds);
    const bool reuseSbfHashes = (seed_bf->get_hash_num_per_seed() == cbf->get_hash_num()) &&
        (seed_bf->get_k() == kMerSize) &&
        (seed_bf->get_seeds() == seeds);

    const std::vector<std::string> singleBatch = {read};
    return flushReadBatch(singleBatch, cbf, seed_bf, parsedCbfSeeds, reuseSbfHashes, kMerSize);
}

void populate(const std::string& seq, const std::vector<std::string>& seeds, btllib::CountingBloomFilter16* cbf)
{
    if (seeds.empty() || cbf == nullptr) return;
    btllib::SeedNtHash h(seq, seeds, cbf->get_hash_num(), seeds[0].size());
    while (h.roll()) {
        for (unsigned i = 0; i < seeds.size(); i++) {
            cbf->insert(h.hashes() + i * cbf->get_hash_num());
        }
    }
}
