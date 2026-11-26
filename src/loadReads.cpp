#include <loadReads.hpp>

namespace {
struct ReadSlice {
	const char* data;
	size_t len;
};

constexpr size_t READ_SLICE_KMER_TARGET = 65536;
constexpr size_t BATCH_MAX_BYTES = 64ULL * 1024ULL * 1024ULL; // 64 MiB max buffered reads per batch

ulong processReadSlice(
	const char* data,
	size_t len,
	CBF *cbf,
	SBF *seed_bf,
	const std::vector<btllib::SpacedSeed>& parsedCbfSeeds,
	unsigned cbfHashNum,
	unsigned sbfHashNum,
	bool reuseSbfHashes,
	uint kMerSize)
{
	ulong hits = 0;
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

void flushReadBatch(
	const std::vector<std::string>& batch,
	CBF *cbf,
	SBF *seed_bf,
	const std::vector<btllib::SpacedSeed>& parsedCbfSeeds,
	bool reuseSbfHashes,
	uint kMerSize)
{
	if (batch.empty()) return;

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

	#pragma omp parallel for schedule(dynamic, 1) shared(slices, cbf, seed_bf, parsedCbfSeeds) reduction(+ : k_reads)
	for (size_t i = 0; i < slices.size(); i++) {
		k_reads += processReadSlice(
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
}
} // namespace

CBF *loadReads(
	string readsFile, SBF *seed_bf, ulong k_num, double fpr, string seed1, uint kMerSize)
{
	// Verify reads file is accessible before allocating the Counting Bloom Filter
	if (readsFile != "-" && readsFile.find("://") == string::npos) {
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

	// Calculate the optimal number of counter slots (m) and convert to bytes for CountingBloomFilter16 (2 bytes per slot)
	size_t cbfSlots = static_cast<size_t>((-1.0) * (static_cast<double>(k_num) * log(fpr)) / (log(2) * log(2)));
	if (cbfSlots < 1) cbfSlots = 1;
	size_t cbfBytes = cbfSlots * sizeof(uint16_t);

	// Calculate and store the optimal number of hash functions (k)
	uint64_t optHashFunctions = static_cast<uint64_t>(round((static_cast<double>(cbfSlots) / k_num) * log(2)));
	if (optHashFunctions < 1) optHashFunctions = 1;

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >CBF Size:    " << cbfSlots << " counters (" << cbfBytes << " bytes)" << endl;

	// Calculate GB for printing: CountingBloomFilter16 uses 2 bytes (uint16_t) per slot
	cout << "\t\t\t[i] >CBF Memory:   " << std::fixed << std::setprecision(4)
		<< static_cast<double>(cbfBytes) / (1024.0 * 1024.0 * 1024.0) << " GB" << endl;

	cout << "\t\t\t[i] >Optimal Hash Functions: " << optHashFunctions << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	// Initialize the CBF with the byte budget and optimized hash count
	CBF* cbf = new CBF(cbfBytes, optHashFunctions);

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

	auto enqueueRead = [&](std::string& seq) {
		if (seq.empty()) return;
		batchBytes += seq.size();
		readBatch.push_back(std::move(seq));
		seq.clear();
		if (readBatch.size() >= batchMaxReads || batchBytes >= BATCH_MAX_BYTES) {
			flushReadBatch(readBatch, cbf, seed_bf, parsedCbfSeeds, reuseSbfHashes, kMerSize);
			readBatch.clear();
			batchBytes = 0;
		}
	};

	// Stream FASTQ or FASTA reads via state-tracked parser (safe for Phred+33 '@'/'+' quality lines and compressed inputs)
	enum class ReadState { HEADER, FASTQ_SEQ, FASTQ_QUAL, FASTA_SEQ };
	ReadState state = ReadState::HEADER;
	string line, content;
	size_t seqLines = 0;
	size_t qualLen = 0;

	while (readLineFromStream(fp, line)) {
		boost::trim_right(line);
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
		}
		else if (state == ReadState::FASTQ_SEQ) {
			if (line[0] == '+') {
				qualLen = 0;
				state = ReadState::FASTQ_QUAL;
			} else {
				content.append(line);
				seqLines++;
			}
		}
		else if (state == ReadState::FASTQ_QUAL) {
			// Consume quality line(s) unconditionally even if Phred+33 starts with '@' (Q=31) or '+' (Q=10)
			qualLen += line.size();
			if (qualLen >= content.size() || seqLines <= 1) {
				enqueueRead(content);
				state = ReadState::HEADER;
			}
		}
		else if (state == ReadState::FASTA_SEQ) {
			if (line[0] == '>') {
				enqueueRead(content);
			} else {
				content.append(line);
			}
		}
	}

	// Enqueue trailing record and flush remaining batch
	if (!content.empty() && (state == ReadState::FASTA_SEQ || state == ReadState::FASTQ_SEQ || state == ReadState::FASTQ_QUAL)) {
		enqueueRead(content);
	}
	if (!readBatch.empty()) {
		flushReadBatch(readBatch, cbf, seed_bf, parsedCbfSeeds, reuseSbfHashes, kMerSize);
		readBatch.clear();
	}

	readsSource.close();

	const double exactCoverage = (k_assembly > 0)
		? (static_cast<double>(k_reads) / static_cast<double>(k_assembly))
		: 0.0;

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >CBF occupancy: " << cbf->get_occupancy() << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	cout << "\t\t\t[i] >----------------------------------------" << endl;
	cout << "\t\t\t[i] >K-MERS IN BLOOM FILTER: " << k_assembly << endl;
	cout << "\t\t\t[i] >K-MERS IN COUNTING BLOOM FILTER: " << k_reads << endl;
	cout << "\t\t\t[i] >BLOOM FILTER COVERAGE DEPTH: " << exactCoverage << endl;
	cout << "\t\t\t[i] >----------------------------------------" << endl;

	return cbf;
}

void populateCBF(
	CBF *cbf, const string& read, SBF *seed_bf, const string& seed1, uint kMerSize)
{
	if (read.size() < kMerSize || cbf == nullptr || seed_bf == nullptr) return;
	const std::vector<std::string> seeds = {seed1};
	const std::vector<btllib::SpacedSeed> parsedCbfSeeds = btllib::parse_seeds(seeds);
	const bool reuseSbfHashes = (seed_bf->get_hash_num_per_seed() == cbf->get_hash_num()) &&
		(seed_bf->get_k() == kMerSize) &&
		(seed_bf->get_seeds() == seeds);

	std::vector<std::string> singleBatch = {read};
	flushReadBatch(singleBatch, cbf, seed_bf, parsedCbfSeeds, reuseSbfHashes, kMerSize);
}

void populate(const std::string& seq, const std::vector<std::string>& seeds, btllib::CountingBloomFilter16 *cbf)
{
	if (seeds.empty() || cbf == nullptr) return;
	btllib::SeedNtHash h(seq, seeds, cbf->get_hash_num(), seeds[0].size());
	while (h.roll()) {
		for (unsigned i = 0; i < seeds.size(); i++) {
			cbf->insert(h.hashes() + i * cbf->get_hash_num());
		}
	}
}
