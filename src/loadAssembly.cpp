#include <loadAssembly.hpp>

bool loadAssembly(const std::string& assemblyFile, SEQLIST& outAssembly)
{
    outAssembly.clear();

    if (assemblyFile != "-" && assemblyFile.find("://") == std::string::npos) {
        std::ifstream checkFile(assemblyFile);
        if (!checkFile.is_open()) {
            std::cerr << "[!] Error: Could not open assembly file: " << assemblyFile << std::endl;
            return false;
        }
    }

    btllib::DataSource assemblySource(assemblyFile);
    FILE* fp = assemblySource;
    if (fp == nullptr) {
        std::cerr << "[!] Error: Could not open assembly stream: " << assemblyFile << std::endl;
        return false;
    }

    std::string line, name, content;

    while (konuseg::readLineFromStream(fp, line)) {
        konuseg::trim(line);
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (!name.empty() && !content.empty()) {
                outAssembly.emplace_back(name, std::move(content));
                content.clear();
            }

            // Strip leading '>' and truncate at first whitespace for valid BED chrom IDs
            std::string header = line.substr(1);
            konuseg::trim(header);
            const auto wsPos = header.find_first_of(" \t");
            name = (wsPos == std::string::npos) ? header : header.substr(0, wsPos);
        } else {
            content.append(line);
        }
    }

    if (!name.empty() && !content.empty()) {
        outAssembly.emplace_back(name, std::move(content));
    }

    assemblySource.close();
    return true;
}
