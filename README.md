# KoNuSeg
[![GitHub Release](https://img.shields.io/github/v/release/BilkentCompGen/KoNuSeg?color=blue&cacheSeconds=60)](https://github.com/BilkentCompGen/KoNuSeg/releases)
[![Last Commit](https://img.shields.io/github/last-commit/BilkentCompGen/KoNuSeg)](https://github.com/BilkentCompGen/KoNuSeg/commits/main)
[![License](https://img.shields.io/github/license/BilkentCompGen/KoNuSeg?cacheSeconds=60)](https://github.com/BilkentCompGen/KoNuSeg/blob/main/LICENSE)
`KoNuSeg` is a C++ tool for fast, efficient repeat detection using k-mer counting and Bloom filters, backed by `btllib`.

---

## Prerequisites & Requirements

Before building `KoNuSeg`, ensure your system meets the following software and library requirements:

* **Compiler**: GCC / G++ version **11** or higher (must support `C++17` and OpenMP).
* **Dependencies**:
  * [btllib](https://github.com/bcgsc/btllib): Library for bioinformatics data structures (Bloom filters, sequence parsing, etc.). This library is included in this repository (version 1.4.9). 
---

## Installation

### Clone the Repository
```bash
git clone git@github.com:BilkentCompGen/KoNuSeg.git
cd KoNuSeg
```

## Compilation & Verification

The project includes a `Makefile` to simplify compilation and testing. Ensure `g++-11` is available at `/usr/bin/g++-11` on your system, or update the `CC` variable in the `Makefile` to match your local compiler path.

### 1. Build the Executable
To compile the optimized production binary:
```bash
make
```

### 2. Sanity check
Before running `KoNuSeg` on your own dataset, execute the included test suite as a sanity check to verify that the binary was built correctly and all bundled dependencies function properly:
```bash
make run_test
```
### 3. To remove generated binary files:
```bash
make clean
```

## Usage

```bash
./bin/konuseg.out -fq <fastq> -fa <fasta> -o <bed> [options]
```

### Required Parameters

| Parameter | Description |
| :--- | :--- |
| `-fq` | Path to input FASTQ file |
| `-fa` | Path to input FASTA file |
| `-o` | Path to output BED file |

### Optional Parameters

| Option | Default | Description |
| :--- | :--- | :--- |
| `-g` | `3000000000` | Genome size |
| `-s` | `10010111001110100101110011101001` | Seed value |
| `-k` | `32` | K-mer length |
| `-fpr` | `0.01` | False positive rate |
| `-w` | `5` | Filtering window size |
| `-minseg` | `500` | Minimum segment length for merging |
| `-c_thresh` | `1` | Count threshold for merging |
| `-anum` | `1500000000` | Distinct number of k-mers in assembly/reference |
| `-knum` | `1500000000` | Distinct number of k-mers in reads |

If `-anum`  and `-knum` are not known, the size of the genome can be used. 

---

### Example Command

```bash
./bin/konuseg.out \
    -fa sample_ref.fa \
    -fq sample_reads.fq \
    -o output_segments.bed \
    -g 100000000 \
    -s 11111111111111111111111111111111 \
    -k 32 \
    -fpr 0.01 \
    -w 5 \
    -minseg 64 \
    -c_thresh 1
```
