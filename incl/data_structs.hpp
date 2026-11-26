
#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <climits>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#if __has_include(<boost/algorithm/string.hpp>)
#include <boost/algorithm/string.hpp>
#else
namespace boost {
inline void trim_right(std::string& s) {
    const auto pos = s.find_last_not_of(" \t\n\r\f\v");
    if (pos == std::string::npos) s.clear();
    else s.erase(pos + 1);
}
inline void trim(std::string& s) {
    trim_right(s);
    const auto pos = s.find_first_not_of(" \t\n\r\f\v");
    if (pos != std::string::npos) s.erase(0, pos);
}
}
#endif
#if __has_include(<boost/container/stable_vector.hpp>)
#include <boost/container/stable_vector.hpp>
template <typename T>
using KoNuSegVector = boost::container::stable_vector<T>;
#else
template <typename T>
using KoNuSegVector = std::vector<T>;
#endif
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>

#include <chrono>


using namespace std;


class Config{

    public:
        string infile_fastq;
        string infile_fasta;
        string outfile_bed;
        
        // Default values (used if the user doesn't provide them in the terminal)
        long long genome_size = 3000000000; 
        string seed1 = "10010111001110100101110011101001";
        uint k_len = 32;
        ulong k_num = 1500000000;
        ulong a_num = 1500000000;
        double fpr = 0.01;
        uint window_size = 5;
        uint seg_len = 500;
        uint count_threshold = 1;


        Config()
        {
            return;
        };

        ~Config()
        {
            return;
        };


        void printConfig()
        {
            cout << "\t\t[i]: Config file" << endl;
            cout << "\t\t[i]:" << infile_fastq << endl;
            cout << "\t\t[i]:" << infile_fasta << endl;
            cout << "\t\t[i]:" << outfile_bed << endl;
            cout << "\t\t[i]:" << genome_size << endl;
            cout << "\t\t[i]:" << seed1 << endl;
            cout << "\t\t[i]:" << k_len << endl;
            cout << "\t\t[i]:" << k_num << endl;
            cout << "\t\t[i]:" << a_num << endl;
            cout << "\t\t[i]:" << fpr << endl;
            cout << "\t\t[i]:" << window_size << endl;
            cout << "\t\t[i]:" << seg_len << endl;
            cout << "\t\t[i]:" << count_threshold << endl;
            cout << "\t\t[i]: //////////////////////" << endl;
        }
};

class Sequence
{

    public:
        string name;
        string data;

        Sequence()
        {
            return;
        };

        Sequence(string n, string d)
        {
            this->name = n;
            this->data = d;
        };

        ~Sequence()
        {
            return;
        };
};

class ContigInfo{
    
    public:
        string contigName;
        uint start, end;

        ContigInfo()
        {
            
        };


        ContigInfo(string n, uint s, uint e)
        {
            this->contigName = n;
            this->start = s;
            this->end = e;
        };

        ~ContigInfo()
        {
            
        };

        
        void print()
        {
            cout << contigName << "= s" << start << "-e" << end << ":" << endl;
        };
};

class ContigRepeatInfo{
    
    public:
        string contigName;
        uint start, end, numRepeats;

        ContigRepeatInfo()
        {
            
        };


        ContigRepeatInfo(string n, uint s, uint e, uint r)
        {
            // cout << "ovlctr" << endl;
            this->contigName = n;
            this->start = s;
            this->end = e;
            this->numRepeats = r;
        };

        ~ContigRepeatInfo()
        {
            
        };

        
        void print()
        {
            cout << contigName << "= s" << start << "-e" << end << ":" << numRepeats << endl;
        };
};


class SegmentInfo {
    public:
    
        std::string segmentName;
        uint start, end;
        float count;
    
        SegmentInfo(){ return;	};
    
        SegmentInfo(string n, uint s, uint e, float c)
        {
            this->segmentName = n;
            this->start = s;
            this->end = e;
            this->count = c;
            return;
        };
    
};
    
    
class WindowInfo {
    
    public:
        float windowDescriptor;

        WindowInfo()
        {
            return;
        };

        WindowInfo(float n)
        {
            this->windowDescriptor = n;
        return;
        };

        ~WindowInfo()
        {
            return;
        };
};




//The function used to get the assembly results in paper
inline float roundCount(float count)
{
    float intPart = floor(count);
    float diff = count - intPart;

    if (count < 0.25){
        return 0.0f;
    }
    else if (count >= 0.25 && count < 0.65){
        return 0.5f;
    }
    else {
        if (diff <= 0.64f){
            return intPart;
        }
        else {
            return intPart + 1.0f;
        }
    }
};


// aliases
typedef KoNuSegVector<Sequence> SEQLIST;
typedef KoNuSegVector<uint> REPLIST;
typedef btllib::SeedBloomFilter SBF;
typedef btllib::CountingBloomFilter16 CBF;
typedef KoNuSegVector<ContigRepeatInfo> CONTIGLIST;
typedef KoNuSegVector<ContigInfo> CONTIGTABLE;
typedef KoNuSegVector<SegmentInfo> SEGMENTINFO;
typedef KoNuSegVector<WindowInfo> WINDOWINFO;
typedef KoNuSegVector< SEGMENTINFO* > SEGMENTTABLE;




inline ulong k_assembly=0;
inline ulong k_reads=0;
inline SEGMENTTABLE* segTable = new SEGMENTTABLE();
