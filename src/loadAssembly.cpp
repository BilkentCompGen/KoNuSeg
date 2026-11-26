#include <loadAssembly.hpp>


SEQLIST *loadAssembly(string assemblyFile, uint k_len)
{
    (void)k_len;

    if (assemblyFile != "-" && assemblyFile.find("://") == string::npos) {
        std::ifstream checkFile(assemblyFile);
        if (!checkFile.is_open()) {
            std::cerr << "[!] Error: Could not open assembly file: " << assemblyFile << std::endl;
            return nullptr;
        }
    }

    btllib::DataSource assemblySource(assemblyFile);
    FILE* fp = assemblySource;
    if (fp == nullptr) {
        std::cerr << "[!] Error: Could not open assembly stream: " << assemblyFile << std::endl;
        return nullptr;
    }

    SEQLIST *myAssembly = new SEQLIST();
    std::string line, name, content;

    while (readLineFromStream(fp, line)) {
        boost::trim(line);
        if (line.empty()) {
            continue;
        }

        // Identify the FASTA header line
        if (line[0] == '>') {
            // If we have DNA accumulated from a previous record, save it now
            if (!name.empty() && !content.empty()) {
                myAssembly->push_back(Sequence(name, content));
                content.clear();
            }

            // Remove the leading '>' and truncate at the first whitespace for valid BED chrom IDs
            std::string header = line.substr(1);
            boost::trim(header);
            const auto wsPos = header.find_first_of(" \t");
            name = (wsPos == std::string::npos) ? header : header.substr(0, wsPos);
        } else {
            content.append(line);
        }
    }

    // Capture the very last sequence in the file
    if (!name.empty() && !content.empty()) {
        myAssembly->push_back(Sequence(name, content));
    }

    assemblySource.close();
    return myAssembly;
}