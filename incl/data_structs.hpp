
#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include <boost/container/list.hpp>
#include <boost/container/stable_vector.hpp>
#include <boost/range/algorithm.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <boost/container/stable_vector.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/statistics/median.hpp>

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
    
        /*
        SegmentInfo(const SegmentInfo& val)
        {
            this->segmentName = val.segmentName;
            this->start = val.start;
            this->end = val.end;
            this->count = val.count;
            return;
        };
    
        ~SegmentInfo() { return; };

        public:
        SegmentInfo& operator=(const SegmentInfo& other) {
            this->segmentName = other.segmentName;
            this->start = other.start;
            this->end = other.end;
            this->count = other.count;
            return *this;
        }
        */
    
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
typedef boost::container::stable_vector<Sequence> SEQLIST;
typedef boost::container::stable_vector<uint> REPLIST;
typedef btllib::SeedBloomFilter SBF;
typedef btllib::CountingBloomFilter16 CBF;
typedef boost::container::stable_vector<ContigRepeatInfo> CONTIGLIST;
typedef boost::container::stable_vector<ContigInfo> CONTIGTABLE;
typedef boost::container::stable_vector<SegmentInfo> SEGMENTINFO;
typedef boost::container::stable_vector<WindowInfo> WINDOWINFO;
typedef boost::container::stable_vector< SEGMENTINFO* > SEGMENTTABLE;




inline ulong k_assembly=0;
inline ulong k_reads=0;
inline SEGMENTTABLE* segTable = new SEGMENTTABLE();
