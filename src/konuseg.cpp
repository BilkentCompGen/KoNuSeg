#include<konuseg.hpp>
#include<limits>
#include<stdexcept>

static bool parseUnsigned64(const string& s, unsigned long long& out, unsigned long long minVal = 0, unsigned long long maxVal = numeric_limits<unsigned long long>::max()) {
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

Config* parseCommandLine(int argc, char **argv) {
    Config *config = new Config(); // Starts with all the default values loaded

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

        // Check for file paths (Required)
        if (arg == "-fq") config->infile_fastq = val;
        else if (arg == "-fa") config->infile_fasta = val;
        else if (arg == "-o") config->outfile_bed = val;

        // Check for optional parameters (Overrides defaults if provided)
        else if (arg == "-g") {
            try {
                size_t pos = 0;
                long long g = stoll(val, &pos);
                if (pos != val.size() || g <= 0) throw invalid_argument("invalid genome size");
                config->genome_size = g;
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
        }
        else if (arg == "-s") {
            config->seed1 = val;
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
        }
        else if (arg == "-anum") {
            if (!parseUnsigned64(val, tmpU64, 1, numeric_limits<ulong>::max())) {
                cerr << "Error: Invalid value for -anum: '" << val << "' (must be a positive 64-bit integer)." << endl;
                delete config;
                return nullptr;
            }
            config->a_num = static_cast<ulong>(tmpU64);
        }
        else {
            cerr << "Error: Unknown option '" << arg << "'." << endl;
            delete config;
            return nullptr;
        }
    }

    // Force user to give the file paths
    if (config->infile_fastq.empty() || config->infile_fasta.empty() || config->outfile_bed.empty()) {
        cerr << "Error: You must provide -fq (fastq), -fa (fasta), and -o (bed) arguments." << endl;
        delete config;
        return nullptr;
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
         << "  -fq               Path to input FASTQ file\n"
         << "  -fa               Path to input FASTA file\n"
         << "  -o                Path to output BED file\n"
         << "Options (with defaults):\n"
         << "  -g                Genome size (default: 3000000000)\n"
         << "  -s                Seed value (default: 10010111001110100101110011101001)\n"
         << "  -k                K-mer length (default: 32)\n"
         << "  -fpr              False positive rate (default: 0.01)\n"
         << "  -w                Filtering window size (default: 5)\n"
         << "  -minseg           Minimum segment length for the merging (default:500)\n"
         << "  -c_thresh         Count threshold for the merging (default:1)\n"
         << "  -anum             Distinct number of kmers in assembly/reference (default: 1500000000)\n" 
         << "  -knum             Distinct number of kmers in reads (default: 1500000000)\n";
}




int main(int argc, char **argv)
{
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds;

	// Check if user provided arguments    
    if (argc <= 1) {
        printUsage(argv[0]); 
        return -1;
    }

    cout << "\t\t[i] > Processing command line arguments..." << endl;
    
    Config *config_file = parseCommandLine(argc, argv);
    if (config_file == nullptr) {
        printUsage(argv[0]);
        return 1;
    }

    config_file->printConfig();

	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////Load Assemly//////////////////////////////////////////////////////////
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Loading assembly from: " << config_file->infile_fasta << endl;
	SEQLIST *assembly = loadAssembly(config_file->infile_fasta, config_file->k_len);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////Populate Bloom Filter////////////////////////////////////////////////
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Loading reads from: " << config_file->infile_fastq << endl;
	cout << "\t\t[i] > Populating Seed Bloom Filter..." << endl;
	SBF *sbf = populateBF(assembly, config_file->seed1, config_file->a_num, config_file->fpr, config_file->k_len);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////Populate Counting Bloom Filter/////////////////////////////////////////
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Loading reads from: " << config_file->infile_fastq << endl;
	cout << "\t\t[i] > Populating Counting Bloom Filter..." << endl;
	CBF *cbf = loadReads(config_file->infile_fastq, config_file->k_len, assembly, sbf, config_file->k_num, config_file->fpr, config_file->seed1, config_file->k_len);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////Querying Repeats////////////////////////////////////////////////////////
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Querying repeats & creating segments..." << endl;
	checkRepeats(assembly, cbf, config_file->k_len, config_file->seed1, config_file->genome_size, config_file->window_size, config_file->outfile_bed, config_file->seg_len, config_file->count_threshold);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	delete sbf;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////Merging Segments////////////////////////////////////////////////////////	
	start = std::chrono::steady_clock::now();
	cout << "\t\t[i] > Merging segments..." << endl;
	mergeSegments(config_file->outfile_bed, config_file->k_len, config_file->seg_len, config_file->count_threshold);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	std::cout << "\t\t[i] > Wall-clock elapsed time: " << elapsed_seconds.count() << "s" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	delete cbf;
	delete assembly;
	delete config_file;
	delete segTable;


	return 0;
}
