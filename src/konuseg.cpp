#include <konuseg.hpp>
#include <limits>
#include <stdexcept>

namespace {
bool parseUnsigned64(const string& s, unsigned long long& out, unsigned long long minVal = 0, unsigned long long maxVal = numeric_limits<unsigned long long>::max()) {
    if (s.empty() || s[0] == '-') return false;
    try {
        size_t pos = 0;
        unsigned long long val = stoull(s, &pos);
        if (pos != s.size() || val < minVal || val > maxVal) return false;
        out = val;
        return true;
    } catch (const exception&) {
        return false;
    }
}

bool saveFilterMeta(const string& filterPath, ulong kAssemblyVal, ulong kReadsVal) {
    ofstream metaOut(filterPath + ".meta");
    if (!metaOut.is_open()) return false;
    metaOut << kAssemblyVal << " " << kReadsVal << "\n";
    return metaOut.good();
}

bool loadFilterMeta(const string& filterPath, ulong& kAssemblyOut, ulong& kReadsOut) {
    ifstream metaIn(filterPath + ".meta");
    if (!metaIn.is_open()) return false;
    return static_cast<bool>(metaIn >> kAssemblyOut >> kReadsOut);
}

void cleanupSegTable() {
    if (segTable != nullptr) {
        for (size_t i = 0; i < segTable->size(); ++i) {
            delete segTable->at(i);
        }
        delete segTable;
        segTable = nullptr;
    }
}
} // namespace

Config* parseCommandLine(int argc, char **argv) {
    Config *config = new Config(); // Starts with all the default values loaded
    bool k_provided = false;
    bool s_provided = false;
    bool g_provided = false;
    bool anum_provided = false;
    bool knum_provided = false;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            delete config;
            return nullptr;
        }

        if (i + 1 >= argc) {
            cerr << "Error: Missing value for argument '" << arg << "'." << endl;
            delete config;
            return nullptr;
        }

        string val = argv[++i];
        unsigned long long tmpU64 = 0;

        // Check for file paths
        if (arg == "-fq") config->infile_fastq = val;
        else if (arg == "-fa") config->infile_fasta = val;
        else if (arg == "-o") config->outfile_bed = val;
        else if (arg == "-save_bf" || arg == "--save-bf") config->save_bf_path = val;
        else if (arg == "-load_bf" || arg == "--load-bf") config->load_bf_path = val;
        else if (arg == "-save_cbf" || arg == "--save-cbf") config->save_cbf_path = val;
        else if (arg == "-load_cbf" || arg == "--load-cbf") config->load_cbf_path = val;

        // Check for optional parameters (Overrides defaults if provided)
        else if (arg == "-g") {
            try {
                size_t pos = 0;
                long long g = stoll(val, &pos);
                if (pos != val.size() || g <= 0) throw invalid_argument("invalid genome size");
                config->genome_size = g;
                g_provided = true;
            } catch (const exception&) {
                cerr << "Error: Invalid value for -g: '" << val << "' (must be a positive integer)." << endl;
                delete config;
                return nullptr;
            }
        }
        else if (arg == "-k") {
            if (!parseUnsigned64(val, tmpU64, 1, numeric_limits<uint>::max())) {
                cerr << "Error: Invalid value for -k: '" << val << "' (must be a positive integer)." << endl;
                delete config;
                return nullptr;
            }
            config->k_len = static_cast<uint>(tmpU64);
            k_provided = true;
        }
        else if (arg == "-s") {
            config->seed1 = val;
            s_provided = true;
        }
        else if (arg == "-fpr") {
            try {
                size_t pos = 0;
                double fpr = stod(val, &pos);
                if (pos != val.size() || !(fpr > 0.0 && fpr < 1.0)) throw invalid_argument("invalid fpr");
                config->fpr = fpr;
            } catch (const exception&) {
                cerr << "Error: Invalid value for -fpr: '" << val << "' (must be in (0, 1))." << endl;
                delete config;
                return nullptr;
            }
        }
        else if (arg == "-w") {
            if (!parseUnsigned64(val, tmpU64, 1, numeric_limits<uint>::max())) {
                cerr << "Error: Invalid value for -w: '" << val << "' (must be >= 1)." << endl;
                delete config;
                return nullptr;
            }
            config->window_size = static_cast<uint>(tmpU64);
        }
        else if (arg == "-minseg") {
            if (!parseUnsigned64(val, tmpU64, 1, numeric_limits<uint>::max())) {
                cerr << "Error: Invalid value for -minseg: '" << val << "' (must be >= 1)." << endl;
                delete config;
                return nullptr;
            }
            config->seg_len = static_cast<uint>(tmpU64);
        }
        else if (arg == "-c_thresh") {
            if (!parseUnsigned64(val, tmpU64, 0, numeric_limits<uint>::max())) {
                cerr << "Error: Invalid value for -c_thresh: '" << val << "' (must be >= 0)." << endl;
                delete config;
                return nullptr;
            }
            config->count_threshold = static_cast<uint>(tmpU64);
        }
        else if (arg == "-knum") {
            if (!parseUnsigned64(val, tmpU64, 1, numeric_limits<ulong>::max())) {
                cerr << "Error: Invalid value for -knum: '" << val << "' (must be a positive 64-bit integer)." << endl;
                delete config;
                return nullptr;
            }
            config->k_num = static_cast<ulong>(tmpU64);
            knum_provided = true;
        }
        else if (arg == "-anum") {
            if (!parseUnsigned64(val, tmpU64, 1, numeric_limits<ulong>::max())) {
                cerr << "Error: Invalid value for -anum: '" << val << "' (must be a positive 64-bit integer)." << endl;
                delete config;
                return nullptr;
            }
            config->a_num = static_cast<ulong>(tmpU64);
            anum_provided = true;
        }
        else {
            cerr << "Error: Unknown option '" << arg << "'." << endl;
            delete config;
            return nullptr;
        }
    }

    // Require -fa and -o always, and require -fq unless loading a prebuilt Counting Bloom Filter (-load_cbf)
    if (config->infile_fasta.empty() || config->outfile_bed.empty() ||
        (config->infile_fastq.empty() && config->load_cbf_path.empty())) {
        cerr << "Error: You must provide -fa (fasta), -o (bed), and either -fq (reads) or -load_cbf (prebuilt CBF)." << endl;
        delete config;
        return nullptr;
    }

    // Use -g (genome_size) as the fallback filter sizing bound when -g is given without -anum / -knum
    if (g_provided && !anum_provided) {
        config->a_num = static_cast<ulong>(config->genome_size);
    }
    if (g_provided && !knum_provided) {
        config->k_num = static_cast<ulong>(config->genome_size);
    }

    // Auto-synchronize -k and -s when only one of the two is explicitly provided
    if (k_provided && !s_provided) {
        config->seed1 = string(config->k_len, '1');
    } else if (s_provided && !k_provided) {
        if (config->seed1.empty() || config->seed1.size() > numeric_limits<uint>::max()) {
            cerr << "Error: Invalid seed length (-s)." << endl;
            delete config;
            return nullptr;
        }
        config->k_len = static_cast<uint>(config->seed1.size());
    }

    // Validate seed string and compatibility with k-mer length
    if (config->seed1.size() != config->k_len) {
        cerr << "Error: Seed length (-s: " << config->seed1.size()
             << ") must match k-mer length (-k: " << config->k_len << ")." << endl;
        delete config;
        return nullptr;
    }
    bool hasCarePos = false;
    for (char c : config->seed1) {
        if (c == '1') hasCarePos = true;
        else if (c != '0') {
            cerr << "Error: Seed (-s) must contain only '0' and '1' characters." << endl;
            delete config;
            return nullptr;
        }
    }
    if (!hasCarePos) {
        cerr << "Error: Seed (-s) must contain at least one '1' position." << endl;
        delete config;
        return nullptr;
    }

    return config;
}


void printUsage(const char* programName) {
    cerr << "Usage: " << programName << " -fq <fastq> -fa <fasta> -o <bed> [options]\n"
         << "Required:\n"
         << "  -fq               Path to input reads (FASTQ/FASTA, supports .gz/.bz2/.xz; optional with -load_cbf)\n"
         << "  -fa               Path to input assembly FASTA file (supports .gz/.bz2/.xz)\n"
         << "  -o                Path to output BED file\n"
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
         << "  -load_cbf         Load prebuilt Counting Bloom Filter from file\n";
}


int main(int argc, char **argv)
{
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds;

	if (argc <= 1) {
		printUsage(argv[0]);
		cleanupSegTable();
		return 1;
	}

	cout << "\t\t[i] > Processing command line arguments..." << endl;

	Config *config_file = parseCommandLine(argc, argv);
	if (config_file == nullptr) {
		printUsage(argv[0]);
		cleanupSegTable();
		return 1;
	}

	config_file->printConfig();

	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;

	/////////////////////////////////////////////////////////////////Load Assembly/////////////////////////////////////////////////////////
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Loading assembly from: " << config_file->infile_fasta << endl;
	SEQLIST *assembly = loadAssembly(config_file->infile_fasta, config_file->k_len);
	if (assembly == nullptr || assembly->empty()) {
		if (assembly != nullptr) {
			cerr << "[!] Error: No valid sequences found in assembly file: " << config_file->infile_fasta << endl;
			delete assembly;
		}
		delete config_file;
		cleanupSegTable();
		return 1;
	}
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;

	SBF *sbf = nullptr;
	CBF *cbf = nullptr;

	//////////////////////////////////////////////////////////Seed & Counting Bloom Filters////////////////////////////////////////////////
	if (!config_file->load_cbf_path.empty()) {
		start = std::chrono::steady_clock::now();
		cout << "\t\t[i] > Loading prebuilt Counting Bloom Filter from: " << config_file->load_cbf_path << endl;
		ifstream testCbf(config_file->load_cbf_path);
		if (!testCbf.is_open() || !CBF::is_bloom_file(config_file->load_cbf_path)) {
			cerr << "[!] Error: Invalid or unreadable Counting Bloom Filter file: " << config_file->load_cbf_path << endl;
			delete assembly;
			delete config_file;
			cleanupSegTable();
			return 1;
		}
		cbf = new CBF(config_file->load_cbf_path);
		if (!loadFilterMeta(config_file->load_cbf_path, k_assembly, k_reads)) {
			cerr << "[!] Warning: Metadata file '" << config_file->load_cbf_path << ".meta' not found; estimating assembly kmers at 1x coverage." << endl;
			k_assembly = countAssemblyKmers(assembly, config_file->k_len);
			k_reads = k_assembly;
		}
		end = std::chrono::steady_clock::now();
		elapsed_seconds = end - start;
		std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	} else {
		// Stage 1: Build or Load Seed Bloom Filter
		start = std::chrono::steady_clock::now();
		if (!config_file->load_bf_path.empty()) {
			cout << "\t\t[i] > Loading prebuilt Seed Bloom Filter from: " << config_file->load_bf_path << endl;
			ifstream testBf(config_file->load_bf_path);
			if (!testBf.is_open() || !SBF::is_bloom_file(config_file->load_bf_path)) {
				cerr << "[!] Error: Invalid or unreadable Seed Bloom Filter file: " << config_file->load_bf_path << endl;
				delete assembly;
				delete config_file;
				cleanupSegTable();
				return 1;
			}
			sbf = new SBF(config_file->load_bf_path);
			config_file->k_len = sbf->get_k();
			if (!sbf->get_seeds().empty()) {
				config_file->seed1 = sbf->get_seeds()[0];
			}
			ulong dummyReads = 0;
			if (!loadFilterMeta(config_file->load_bf_path, k_assembly, dummyReads)) {
				k_assembly = countAssemblyKmers(assembly, config_file->k_len);
			}
		} else {
			cout << "\t\t[i] > Populating Seed Bloom Filter..." << endl;
			sbf = populateBF(assembly, config_file->seed1, config_file->a_num, config_file->fpr, config_file->k_len);
			if (sbf == nullptr) {
				delete assembly;
				delete config_file;
				cleanupSegTable();
				return 1;
			}
			if (!config_file->save_bf_path.empty()) {
				cout << "\t\t[i] > Saving Seed Bloom Filter to: " << config_file->save_bf_path << endl;
				sbf->save(config_file->save_bf_path);
				saveFilterMeta(config_file->save_bf_path, k_assembly, 0);
			}
		}
		end = std::chrono::steady_clock::now();
		elapsed_seconds = end - start;
		std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;

		// Stage 2: Stream Reads & Populate Counting Bloom Filter
		start = std::chrono::steady_clock::now();
		cout << "\t\t[i] > Loading reads from: " << config_file->infile_fastq << endl;
		cout << "\t\t[i] > Populating Counting Bloom Filter..." << endl;
		cbf = loadReads(config_file->infile_fastq, sbf, config_file->k_num, config_file->fpr, config_file->seed1, config_file->k_len);
		if (cbf == nullptr) {
			delete sbf;
			delete assembly;
			delete config_file;
			cleanupSegTable();
			return 1;
		}
		if (!config_file->save_cbf_path.empty()) {
			cout << "\t\t[i] > Saving Counting Bloom Filter to: " << config_file->save_cbf_path << endl;
			cbf->save(config_file->save_cbf_path);
			saveFilterMeta(config_file->save_cbf_path, k_assembly, k_reads);
		}
		end = std::chrono::steady_clock::now();
		elapsed_seconds = end - start;
		std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	}

	delete sbf;
	sbf = nullptr;

	////////////////////////////////////////////////////////////////Querying Repeats////////////////////////////////////////////////////////
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Querying repeats & creating segments..." << endl;
	bool okRepeats = checkRepeats(assembly, cbf, config_file->k_len, config_file->seed1, config_file->window_size, config_file->outfile_bed, config_file->seg_len, config_file->count_threshold);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	if (!okRepeats) {
		delete cbf;
		delete assembly;
		delete config_file;
		cleanupSegTable();
		return 1;
	}

	///////////////////////////////////////////////////////////////Merging Segments////////////////////////////////////////////////////////	
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Merging segments..." << endl;
	bool okMerge = mergeSegments(config_file->outfile_bed, config_file->k_len, config_file->seg_len, config_file->count_threshold);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;

	delete cbf;
	delete assembly;
	delete config_file;
	cleanupSegTable();

	return okMerge ? 0 : 1;
}
