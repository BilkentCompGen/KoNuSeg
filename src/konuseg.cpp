#include<konuseg.hpp>

Config* parseCommandLine(int argc, char **argv) {
    Config *config = new Config(); // Starts with all the default values loaded
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        
        // Check for file paths (Required)
        if (arg == "-fq" && i + 1 < argc) config->infile_fastq = argv[++i];
        else if (arg == "-fa" && i + 1 < argc) config->infile_fasta = argv[++i];
        else if (arg == "-o" && i + 1 < argc) config->outfile_bed = argv[++i];
        
        // Check for optional parameters (Overrides defaults if provided)
        else if (arg == "-g" && i + 1 < argc) config->genome_size = stoll(argv[++i]);
        else if (arg == "-k" && i + 1 < argc) config->k_len = stoi(argv[++i]);
		else if (arg == "-s" && i + 1 < argc) config->seed1 = argv[++i];
        else if (arg == "-fpr" && i + 1 < argc) config->fpr = stod(argv[++i]);
        else if (arg == "-w" && i + 1 < argc) config->window_size = stoi(argv[++i]);
        else if (arg == "-minseg" && i + 1 < argc) config->seg_len = stoi(argv[++i]);
        else if (arg == "-c_thresh" && i + 1 < argc) config->count_threshold = stoi(argv[++i]);
        else if (arg == "-knum" && i + 1 < argc) config->k_num = stoi(argv[++i]);
        else if (arg == "-anum" && i + 1 < argc) config->a_num = stoi(argv[++i]);
    }

    // Force user to give the file paths
    if (config->infile_fastq.empty() || config->infile_fasta.empty() || config->outfile_bed.empty()) {
        cerr << "Error: You must provide -fq (fastq), -fa (fasta), and -o (bed) arguments." << endl;
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
    if (config_file == nullptr) return -1;

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
