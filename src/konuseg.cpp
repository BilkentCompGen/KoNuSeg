#include <konuseg.hpp>
#include <limits>
#include <stdexcept>

namespace {

enum class ParseResult {
    OK,
    HELP,
    ERROR
};

bool parseUnsigned64(
    const std::string& s,
    unsigned long long& out,
    unsigned long long minVal = 0,
    unsigned long long maxVal = std::numeric_limits<unsigned long long>::max())
{
    if (s.empty() || s[0] == '-') return false;
    try {
        size_t pos = 0;
        const unsigned long long val = std::stoull(s, &pos);
        if (pos != s.size() || val < minVal || val > maxVal) return false;
        out = val;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool saveFilterMeta(const std::string& filterPath, uint64_t kAssemblyVal, uint64_t kReadsVal) {
    std::ofstream metaOut(filterPath + ".meta");
    if (!metaOut.is_open()) return false;
    metaOut << kAssemblyVal << " " << kReadsVal << "\n";
    return metaOut.good();
}

bool loadFilterMeta(const std::string& filterPath, uint64_t& kAssemblyOut, uint64_t& kReadsOut) {
    std::ifstream metaIn(filterPath + ".meta");
    if (!metaIn.is_open()) return false;
    return static_cast<bool>(metaIn >> kAssemblyOut >> kReadsOut);
}

} // namespace

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " -fq <fastq> -fa <fasta> -o <bed> [options]\n"
              << "Required:\n"
              << "  -fq               Path to input reads (FASTQ/FASTA, supports .gz/.bz2/.xz; optional with -load_cbf)\n"
              << "  -fa               Path to input assembly FASTA file (supports .gz/.bz2/.xz)\n"
              << "  -o                Path to output BED file (merged output written to <stem>_merged.bed)\n"
              << "Options (with defaults):\n"
              << "  -g                Genome size; also sets -anum/-knum when omitted (default: 3000000000)\n"
              << "  -s                Spaced seed mask (default: 10010111001110100101110011101001)\n"
              << "  -k                K-mer length (default: 32)\n"
              << "  -fpr              Target Bloom filter false-positive rate (default: 0.01)\n"
              << "  -w                Sliding window size for median outlier filtering (default: 5)\n"
              << "  -minseg           Minimum segment length for Stage 2 merging (default: 500)\n"
              << "  -c_thresh         Copy-number threshold for segmentation & merging compatibility (default: 1)\n"
              << "  -anum             Expected distinct kmers in assembly (default: 1500000000, or -g if given)\n"
              << "  -knum             Expected distinct kmers in reads (default: 1500000000, or -g if given)\n"
              << "  -save_bf          Save populated Seed Bloom Filter to file\n"
              << "  -load_bf          Load prebuilt Seed Bloom Filter from file\n"
              << "  -save_cbf         Save populated Counting Bloom Filter to file\n"
              << "  -load_cbf         Load prebuilt Counting Bloom Filter from file\n"
              << "  -h, --help        Display this help message and exit\n";
}

ParseResult parseCommandLine(int argc, char** argv, Config& config) {
    bool k_provided = false;
    bool s_provided = false;
    bool g_provided = false;
    bool anum_provided = false;
    bool knum_provided = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            return ParseResult::HELP;
        }

        if (i + 1 >= argc) {
            std::cerr << "Error: Missing value for argument '" << arg << "'." << std::endl;
            return ParseResult::ERROR;
        }

        const std::string val = argv[++i];
        unsigned long long tmpU64 = 0;

        if (arg == "-fq") config.infile_fastq = val;
        else if (arg == "-fa") config.infile_fasta = val;
        else if (arg == "-o") config.outfile_bed = val;
        else if (arg == "-save_bf" || arg == "--save-bf") config.save_bf_path = val;
        else if (arg == "-load_bf" || arg == "--load-bf") config.load_bf_path = val;
        else if (arg == "-save_cbf" || arg == "--save-cbf") config.save_cbf_path = val;
        else if (arg == "-load_cbf" || arg == "--load-cbf") config.load_cbf_path = val;
        else if (arg == "-g") {
            try {
                size_t pos = 0;
                const long long g = std::stoll(val, &pos);
                if (pos != val.size() || g <= 0) throw std::invalid_argument("invalid genome size");
                config.genome_size = g;
                g_provided = true;
            } catch (const std::exception&) {
                std::cerr << "Error: Invalid value for -g: '" << val << "' (must be a positive integer)." << std::endl;
                return ParseResult::ERROR;
            }
        }
        else if (arg == "-k") {
            if (!parseUnsigned64(val, tmpU64, 1, std::numeric_limits<uint32_t>::max())) {
                std::cerr << "Error: Invalid value for -k: '" << val << "' (must be a positive integer)." << std::endl;
                return ParseResult::ERROR;
            }
            config.k_len = static_cast<uint32_t>(tmpU64);
            k_provided = true;
        }
        else if (arg == "-s") {
            config.seed1 = val;
            s_provided = true;
        }
        else if (arg == "-fpr") {
            try {
                size_t pos = 0;
                const double fpr = std::stod(val, &pos);
                if (pos != val.size() || !(fpr > 0.0 && fpr < 1.0)) throw std::invalid_argument("invalid fpr");
                config.fpr = fpr;
            } catch (const std::exception&) {
                std::cerr << "Error: Invalid value for -fpr: '" << val << "' (must be in (0, 1))." << std::endl;
                return ParseResult::ERROR;
            }
        }
        else if (arg == "-w") {
            if (!parseUnsigned64(val, tmpU64, 1, std::numeric_limits<uint32_t>::max())) {
                std::cerr << "Error: Invalid value for -w: '" << val << "' (must be >= 1)." << std::endl;
                return ParseResult::ERROR;
            }
            config.window_size = static_cast<uint32_t>(tmpU64);
        }
        else if (arg == "-minseg") {
            if (!parseUnsigned64(val, tmpU64, 1, std::numeric_limits<uint32_t>::max())) {
                std::cerr << "Error: Invalid value for -minseg: '" << val << "' (must be >= 1)." << std::endl;
                return ParseResult::ERROR;
            }
            config.seg_len = static_cast<uint32_t>(tmpU64);
        }
        else if (arg == "-c_thresh") {
            if (!parseUnsigned64(val, tmpU64, 0, std::numeric_limits<uint32_t>::max())) {
                std::cerr << "Error: Invalid value for -c_thresh: '" << val << "' (must be >= 0)." << std::endl;
                return ParseResult::ERROR;
            }
            config.count_threshold = static_cast<uint32_t>(tmpU64);
        }
        else if (arg == "-knum") {
            if (!parseUnsigned64(val, tmpU64, 1, std::numeric_limits<uint64_t>::max())) {
                std::cerr << "Error: Invalid value for -knum: '" << val << "' (must be a positive 64-bit integer)." << std::endl;
                return ParseResult::ERROR;
            }
            config.k_num = static_cast<uint64_t>(tmpU64);
            knum_provided = true;
        }
        else if (arg == "-anum") {
            if (!parseUnsigned64(val, tmpU64, 1, std::numeric_limits<uint64_t>::max())) {
                std::cerr << "Error: Invalid value for -anum: '" << val << "' (must be a positive 64-bit integer)." << std::endl;
                return ParseResult::ERROR;
            }
            config.a_num = static_cast<uint64_t>(tmpU64);
            anum_provided = true;
        }
        else {
            std::cerr << "Error: Unknown option '" << arg << "'." << std::endl;
            return ParseResult::ERROR;
        }
    }

    if (config.infile_fasta.empty() || config.outfile_bed.empty() ||
        (config.infile_fastq.empty() && config.load_cbf_path.empty())) {
        std::cerr << "Error: You must provide -fa (fasta), -o (bed), and either -fq (reads) or -load_cbf (prebuilt CBF)." << std::endl;
        return ParseResult::ERROR;
    }

    if (g_provided && !anum_provided) {
        config.a_num = static_cast<uint64_t>(config.genome_size);
    }
    if (g_provided && !knum_provided) {
        config.k_num = static_cast<uint64_t>(config.genome_size);
    }

    if (k_provided && !s_provided) {
        config.seed1 = std::string(config.k_len, '1');
    } else if (s_provided && !k_provided) {
        if (config.seed1.empty() || config.seed1.size() > std::numeric_limits<uint32_t>::max()) {
            std::cerr << "Error: Invalid seed length (-s)." << std::endl;
            return ParseResult::ERROR;
        }
        config.k_len = static_cast<uint32_t>(config.seed1.size());
    }

    if (config.seed1.size() != config.k_len) {
        std::cerr << "Error: Seed length (-s: " << config.seed1.size()
                  << ") must match k-mer length (-k: " << config.k_len << ")." << std::endl;
        return ParseResult::ERROR;
    }

    bool hasCarePos = false;
    for (char c : config.seed1) {
        if (c == '1') hasCarePos = true;
        else if (c != '0') {
            std::cerr << "Error: Seed (-s) must contain only '0' and '1' characters." << std::endl;
            return ParseResult::ERROR;
        }
    }
    if (!hasCarePos) {
        std::cerr << "Error: Seed (-s) must contain at least one '1' position." << std::endl;
        return ParseResult::ERROR;
    }

    return ParseResult::OK;
}

int main(int argc, char** argv)
{
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds{};

    if (argc <= 1) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    Config config;
    const ParseResult parseStatus = parseCommandLine(argc, argv, config);
    if (parseStatus == ParseResult::HELP) {
        printUsage(argv[0]);
        return EXIT_SUCCESS;
    }
    if (parseStatus == ParseResult::ERROR) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    std::cout << "\t\t[i] > Processing command line arguments..." << std::endl;
    config.printConfig();

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

    /////////////////////////////////////////////////////////////////Load Assembly/////////////////////////////////////////////////////////
    start = std::chrono::steady_clock::now();
    std::cout << "\t\t[i] > Loading assembly from: " << config.infile_fasta << std::endl;
    SEQLIST assembly;
    if (!loadAssembly(config.infile_fasta, assembly)) {
        return EXIT_FAILURE;
    }
    if (assembly.empty()) {
        std::cerr << "[!] Error: No valid sequences found in assembly file: " << config.infile_fasta << std::endl;
        return EXIT_FAILURE;
    }
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

    PipelineState state;
    std::unique_ptr<SBF> sbf;
    std::unique_ptr<CBF> cbf;

    //////////////////////////////////////////////////////////Seed & Counting Bloom Filters////////////////////////////////////////////////
    if (!config.load_cbf_path.empty()) {
        start = std::chrono::steady_clock::now();
        std::cout << "\t\t[i] > Loading prebuilt Counting Bloom Filter from: " << config.load_cbf_path << std::endl;
        std::ifstream testCbf(config.load_cbf_path);
        if (!testCbf.is_open() || !CBF::is_bloom_file(config.load_cbf_path)) {
            std::cerr << "[!] Error: Invalid or unreadable Counting Bloom Filter file: " << config.load_cbf_path << std::endl;
            return EXIT_FAILURE;
        }
        cbf = std::make_unique<CBF>(config.load_cbf_path);
        if (!loadFilterMeta(config.load_cbf_path, state.k_assembly, state.k_reads)) {
            std::cerr << "[!] Warning: Metadata file '" << config.load_cbf_path
                      << ".meta' not found; estimating assembly kmers at 1x coverage." << std::endl;
            state.k_assembly = countAssemblyKmers(assembly, config.k_len);
            state.k_reads = state.k_assembly;
        }
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end - start;
        std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
    } else {
        start = std::chrono::steady_clock::now();
        if (!config.load_bf_path.empty()) {
            std::cout << "\t\t[i] > Loading prebuilt Seed Bloom Filter from: " << config.load_bf_path << std::endl;
            std::ifstream testBf(config.load_bf_path);
            if (!testBf.is_open() || !SBF::is_bloom_file(config.load_bf_path)) {
                std::cerr << "[!] Error: Invalid or unreadable Seed Bloom Filter file: " << config.load_bf_path << std::endl;
                return EXIT_FAILURE;
            }
            sbf = std::make_unique<SBF>(config.load_bf_path);
            config.k_len = sbf->get_k();
            if (!sbf->get_seeds().empty()) {
                config.seed1 = sbf->get_seeds()[0];
            }
            uint64_t dummyReads = 0;
            if (!loadFilterMeta(config.load_bf_path, state.k_assembly, dummyReads)) {
                state.k_assembly = countAssemblyKmers(assembly, config.k_len);
            }
        } else {
            std::cout << "\t\t[i] > Populating Seed Bloom Filter..." << std::endl;
            sbf = populateBF(assembly, config.seed1, config.a_num, config.fpr, config.k_len, state.k_assembly);
            if (!sbf) {
                return EXIT_FAILURE;
            }
            if (!config.save_bf_path.empty()) {
                std::cout << "\t\t[i] > Saving Seed Bloom Filter to: " << config.save_bf_path << std::endl;
                sbf->save(config.save_bf_path);
                saveFilterMeta(config.save_bf_path, state.k_assembly, 0);
            }
        }
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end - start;
        std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

        start = std::chrono::steady_clock::now();
        std::cout << "\t\t[i] > Loading reads from: " << config.infile_fastq << std::endl;
        std::cout << "\t\t[i] > Populating Counting Bloom Filter..." << std::endl;
        cbf = loadReads(config.infile_fastq, sbf.get(), config.k_num, config.fpr, config.seed1, config.k_len, state.k_assembly, state.k_reads);
        if (!cbf) {
            return EXIT_FAILURE;
        }
        if (!config.save_cbf_path.empty()) {
            std::cout << "\t\t[i] > Saving Counting Bloom Filter to: " << config.save_cbf_path << std::endl;
            cbf->save(config.save_cbf_path);
            saveFilterMeta(config.save_cbf_path, state.k_assembly, state.k_reads);
        }
        end = std::chrono::steady_clock::now();
        elapsed_seconds = end - start;
        std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
    }

    sbf.reset();

    ////////////////////////////////////////////////////////////////Querying Repeats////////////////////////////////////////////////////////
    start = std::chrono::steady_clock::now();
    std::cout << "\t\t[i] > Querying repeats & creating segments..." << std::endl;
    const bool okRepeats = checkRepeats(
        assembly, cbf.get(), config.k_len, config.seed1,
        config.window_size, config.outfile_bed, config.seg_len, config.count_threshold, state);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
    if (!okRepeats) {
        return EXIT_FAILURE;
    }

    ///////////////////////////////////////////////////////////////Merging Segments////////////////////////////////////////////////////////
    start = std::chrono::steady_clock::now();
    std::cout << "\t\t[i] > Merging segments..." << std::endl;
    const bool okMerge = mergeSegments(config.outfile_bed, config.k_len, config.seg_len, config.count_threshold, state.seg_table);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

    return okMerge ? EXIT_SUCCESS : EXIT_FAILURE;
}
