#include <loadAssembly.hpp>


SEQLIST *loadAssembly(string assemblyFile, uint k_len)
{

	std::string line, name, content;
    std::ifstream assembly(assemblyFile);

    if (!assembly.is_open()) {
        std::cerr << "[!] Error: Could not open assembly file: " << assemblyFile << std::endl;
        return nullptr; 
    }

    SEQLIST *myAssembly = new SEQLIST();

    while (std::getline(assembly, line)) {
        // Handle empty lines or line endings (\r)
        if (line.empty() || line[0] == '\r') {
            continue;
        }

        // Identify the header line
        if (line[0] == '>') {
            // If we have DNA accumulated from a previous record, save it now
            if (!content.empty()) {
                myAssembly->push_back(Sequence(name, content));
                content.clear(); // Reset for the next sequence
            }
            
            // Remove the '>' character
            line.erase(0, 0); 

            // Trim both ends to remove \r, \n, or accidental spaces
            boost::trim(line); 
            name = line;

        } else {
            // 3. This is a sequence line
            boost::trim(line);
            content.append(line);
        }
    }

    // Capture the very last sequence in the file
    if (!content.empty()) {
        myAssembly->push_back(Sequence(name, content));
    }

    assembly.close();
    
    // Return the pointer to the populated list
    return myAssembly;
}