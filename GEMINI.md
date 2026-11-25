# KoNuSeg — Project Knowledgebase

## 1. Overview & Scientific Purpose

**KoNuSeg** is a C++17 / OpenMP bioinformatics tool developed by the Alkan Lab (Bilkent University) and Birol Lab (BC Cancer Research Institute / UBC) for fast, space-efficient, mapping-free detection of repetitive sequences and collapsed segmental duplications in *de novo* genome assemblies.

- **Reference Manuscript**: Zambaku K., Roman-Brenes R., Coombe L., Warren R. L., Wong J., Birol İ., & Alkan C. (2026). *"Space-efficient detection of repeats in genome assemblies using Bloom filters."*
- **Core Idea**: Standard *de novo* repeat annotators inspect only the assembled contigs and therefore cannot distinguish a **collapsed duplication** (assembled as a single copy) from a true single-copy locus. Sequencing reads sample all genomic copies proportionally, so k-mers originating from collapsed repeats appear at multiples of the baseline haploid coverage. Instead of expensive multi-mapping read alignments, KoNuSeg uses a **spaced-seed Bloom filter** (`btllib::SeedBloomFilter`) and a **16-bit Counting Bloom filter** (`btllib::CountingBloomFilter16`) to count error-filtered read k-mers along assembly contigs, smooth hash-collision outliers, and partition contigs into copy-number-annotated BED intervals.
- **Primary Downstream Use Cases**:
  1. **Repeat-aware scaffolding**: Masking segments with estimated copy number exceeding a chosen threshold (`>1`, `>2`, `>5`, or `>10`) via `bedtools maskfasta` before scaffolding draft contigs (e.g., GoldRush contigs scaffolded with `ntLink`) to prevent repeat-induced misjoins.
  2. **Self-alignment repeat profiling**: Supplying an assembly FASTA as both contig and read input (converted to FASTQ) to detect exact/near-exact intragenomic duplications at ~1× baseline depth without whole-genome BLAST alignments.

---

## 2. Repository Structure & External API Policy

> **CRITICAL RULE — DO NOT MODIFY `btllib`**:
> `btllib` (`incl/btllib/`, `incl/cpptoml/`, `incl/tests/`, `lib/libbtllib.a`, and `lib/btllib/`) is a vendored external API (v1.4.9) and **MUST NOT be modified under any circumstances**. All bug fixes, optimizations, and feature changes must be implemented strictly in KoNuSeg's own source and header files (`src/*.cpp`, `incl/*.hpp` outside `incl/btllib/`, and `Makefile`).

```text
KoNuSeg/
├── Makefile                  # Build rules (build, debug, clean, run_test)
├── README.md                 # User installation and CLI reference
├── LICENSE                   # BSD 3-Clause License
├── bin/
│   └── konuseg.out           # Compiled executable target
├── incl/
│   ├── konuseg.hpp           # Umbrella header included by konuseg.cpp
│   ├── data_structs.hpp      # Config, Sequence, SegmentInfo, typedefs, globals, roundCount()
│   ├── loadAssembly.hpp      # FASTA loader declaration
│   ├── populateBF.hpp        # SeedBloomFilter builder declaration
│   ├── loadReads.hpp         # FASTQ streamer & CountingBloomFilter16 builder declarations
│   ├── checkRepeats.hpp      # Contig querying, outlier filtering, Stage 1 segmentation
│   ├── mergeSegments.hpp     # Stage 2 segment pooling & copy-number estimation
│   ├── btllib/               # Bundled btllib v1.4.9 headers (Bloom filters, ntHash, SeqReader, sdsl)
│   ├── cpptoml/              # Header-only TOML parser used by btllib serialization
│   └── tests/                # Bundled upstream btllib test files
├── lib/
│   ├── libbtllib.a           # Prebuilt static btllib v1.4.9 archive
│   └── btllib/               # Upstream btllib Python wrapper artifacts
├── src/
│   ├── konuseg.cpp           # CLI parser, pipeline driver, timing logs, main()
│   ├── loadAssembly.cpp      # Loads contigs from FASTA into in-memory SEQLIST
│   ├── populateBF.cpp        # Populates SeedBloomFilter from assembly k-mers
│   ├── loadReads.cpp         # Streams FASTQ reads and populates CountingBloomFilter16
│   ├── checkRepeats.cpp      # Queries CBF, runs noiseFilter() and createSegments()
│   └── mergeSegments.cpp     # Pools short segments into *_merged.bed via processData()
└── test/
    ├── test.fa               # Small synthetic 5-contig FASTA fixture (seq1: 1x, seq2/seq3: 2x)
    └── test.fq               # Matching 5-read FASTQ fixture for sanity checking
```

---

## 3. Pipeline Architecture & Data Flow

The pipeline executes five sequential phases orchestrated by `main()` in `src/konuseg.cpp`:

### Stage 0: Configuration & Assembly Loading (`src/konuseg.cpp`, `src/loadAssembly.cpp`)
- `parseCommandLine()` populates a `Config` instance (`incl/data_structs.hpp`).
- `loadAssembly()` reads all FASTA records into `SEQLIST` (`boost::container::stable_vector<Sequence>`) in RAM so contigs can be reused during both Bloom filter population and repeat querying.

### Stage 1: Populating the Seed Bloom Filter (`src/populateBF.cpp`)
- Computes optimal filter parameters from expected distinct assembly k-mers (`-anum`, $n$) and false-positive rate (`-fpr`, $p$):
  $$m = -\frac{n \ln p}{(\ln 2)^2}, \qquad \text{bytes} = \left\lceil \frac{m}{8} \right\rceil, \qquad k_{\text{hash}} = \max\!\left(1, \operatorname{round}\!\left(\frac{m}{n} \ln 2\right)\right)$$
- Initializes `btllib::SeedBloomFilter` (`SBF`) with `bfBytes` and the user's spaced seed (`-s`).
- Extracts all `k`-mers across each contig using fine-grained OpenMP parallelism (`#pragma omp parallel for` over k-mer start positions `j`), skipping `N`/non-`ACGT` k-mers via `btllib::SeedNtHash::roll()`, inserting valid spaced-seed hashes into `seed_bf`, and accumulating valid assembly k-mers in global `k_assembly`.

### Stage 2: Streaming Reads & Populating the Counting Bloom Filter (`src/loadReads.cpp`)
- Computes optimal counter slot count `cbfSlots` ($m$) and `optHashFunctions` from expected distinct read k-mers (`-knum`) and `-fpr`, converting to `cbfBytes = cbfSlots * sizeof(uint16_t)` to initialize `btllib::CountingBloomFilter16` (`CBF`, 16-bit atomic counters saturating at `65,535`).
- Streams the input FASTQ file record-by-record without retaining raw reads in RAM.
- For each read (`populateCBF()`), threads extract k-mers in parallel and query `seed_bf->contains(kmer)`:
  - **Invalid (`N`) or Absent in SBF**: Guarded via `!hit_seeds.empty()`; k-mers with `N`s or sequencing errors (absent from the assembly) are discarded.
  - **Present in SBF**: Hashed with `btllib::SeedNtHash` using the same spaced seed (`populate()`) and inserted into `cbf`, incrementing global `k_reads`.
- **Bloom Filter Coverage** is defined as `k_reads / k_assembly` (the effective k-mer depth of coverage in the Counting Bloom filter).

### Stage 3: Querying, Outlier Filtering, & Unmerged Segmentation (`src/checkRepeats.cpp`)
- Iterates over contigs in parallel (`#pragma omp parallel for` over contigs in `checkRepeats()`).
- **Querying (`query()`)**: Re-extracts each contig k-mer, hashes with `SeedNtHash` (returning `0` for `N`-containing k-mers), queries `cbf->contains()`, and stores non-zero k-mer intervals in a thread-local `CONTIGLIST`.
- **Outlier Filtering (`noiseFilter()`)**: Slides a centered window of size `-w` (`windowSize`) out-of-place over the contig's raw k-mer counts (`std::vector<uint>`) to replace hash-collision spikes with the local window median.
- **Unmerged Segmentation (`createSegments()`)**: Inside an `#pragma omp critical` section:
  - Normalizes k-mer counts by continuous floating-point `coverage_depth = float(k_reads) / float(k_assembly)` (guarded when `k_assembly == 0` or `k_reads == 0`).
  - Iterates across all contig k-mers (`i = 1 .. size - 1`), extending the current segment while consecutive k-mers remain within `±count_threshold` (`-c_thresh`) of the segment's running mean copy number (`wMedian / kCount`), have copy number `> 0.1`, and gap `<= k_len`.
  - Writes finalized high-resolution ("unmerged") segments (including the trailing segment at the end of each contig) to `<outfile_bed>` and appends per-contig `SEGMENTINFO*` vectors to global `segTable`.
- `seed_bf` is freed immediately after `checkRepeats()` returns (`src/konuseg.cpp:118`).

### Stage 4: Segment Pooling & Merging (`src/mergeSegments.cpp`)
- `mergeSegments()` iterates sequentially through `segTable` to smooth short segments into longer genomic blocks:
  - Concatenates consecutive unmerged segments on the same contig while accumulated segment length `< min_seg_len` (`-minseg`) and inter-segment gap `<= 3 * k_len`.
  - Once `min_seg_len` is reached (or a contig/gap boundary occurs), `processData()` sorts the pooled segment copy numbers, removes the highest 25% (`size / 4`) to suppress residual hash-collision inflation, computes the median of the remaining 75%, passes it through `roundCount()`, and writes the merged intervals to `<outfile_bed>_merged.bed`.

### Copy-Number Rounding Heuristic (`roundCount()` in `incl/data_structs.hpp`)
- `< 0.25` -> `0.0`
- `[0.25, 0.65)` -> `0.5`
- `>= 0.65`: rounds down to `floor(x)` if fractional part `<= 0.64`, otherwise rounds up to `floor(x) + 1.0`.
- Segments with representative haploid copy state `<= 1` are interpreted as single-copy/neutral; segments `> 1` are classified as duplicated.

---

## 4. Command-Line Interface & Parameter Guidelines

```bash
./bin/konuseg.out -fq <reads.fastq> -fa <contigs.fasta> -o <output.bed> [options]
```

| Flag | Required | Default | Internal Field | Description |
| :--- | :---: | :--- | :--- | :--- |
| `-fq` | Yes | — | `Config::infile_fastq` | Input sequencing reads in uncompressed FASTQ format |
| `-fa` | Yes | — | `Config::infile_fasta` | Input draft assembly contigs in uncompressed FASTA format |
| `-o` | Yes | — | `Config::outfile_bed` | Output path for unmerged BED (merged written to `<o>_merged.bed`) |
| `-g` | No | `3000000000` | `Config::genome_size` | Estimated genome size (currently unused in calculations) |
| `-s` | No | `10010111001110100101110011101001` | `Config::seed1` | Binary spaced-seed mask (`1` = match, `0` = wildcard; length must match `-k`) |
| `-k` | No | `32` | `Config::k_len` | K-mer length (must equal `seed1.size()`) |
| `-fpr` | No | `0.01` | `Config::fpr` | Target Bloom filter false-positive rate |
| `-w` | No | `5` | `Config::window_size` | Sliding window size for median outlier filtering |
| `-minseg` | No | `500` | `Config::seg_len` | Minimum segment length (bp) target during Stage 2 merging |
| `-c_thresh` | No | `1` | `Config::count_threshold` | Copy-number deviation threshold during Stage 1 unmerged segmentation |
| `-anum` | No | `1500000000` | `Config::a_num` | Expected distinct k-mers in the assembly (`SeedBloomFilter` sizing) |
| `-knum` | No | `1500000000` | `Config::k_num` | Expected distinct k-mers in the reads (`CountingBloomFilter16` sizing) |

### Practical Parameter Guidance (from Manuscript Benchmarks)
- **Compact / moderately repetitive genomes** (*C. elegans* ~100 Mbp, *D. melanogaster* ~140 Mbp): Permissive-to-moderate masking thresholds (masking unmerged segments with copy estimate `>1`, `>2`, or `>5`) substantially improve scaffolding NGA50 (up to 4.72× in *C. elegans* and 2.04× in *D. melanogaster*) while reducing misassemblies.
- **Large / repeat-rich genomes** (*S. lycopersicum* ~830 Mbp / 60% repeats, *H. sapiens* ~3.0 Gbp / 51% repeats): Permissive masking (`>1`) over-fragments long repeat paths; use strict copy thresholds (`>5` or `>10`) to reduce misassemblies while preserving baseline contiguity.
- **Thread Scaling**: Parallelism is controlled via standard OpenMP environment variables (e.g., `OMP_NUM_THREADS=32 ./bin/konuseg.out ...`).

---

## 5. Build, Dependencies, & Verification

### Dependencies
1. **GCC 11 (`g++-11`) with C++17 and OpenMP**: **Specifically GCC 11 (`g++-11`) is required to link against the vendored `lib/libbtllib.a`**, because `lib/libbtllib.a` was compiled as slim LTO bytecode (`.gnu.lto_*`, LTO version 11.3, empty `.text` sections). Newer GCC versions (e.g., GCC 13 with LTO 13.1 on Ubuntu 24.04) will compile all `src/*.cpp` translation units cleanly (`g++ -c`), but the linker (`lto1`) will reject GCC 11.3 LTO bytecode (`bytecode stream in file './lib/libbtllib.a' generated with LTO version 11.3 instead of the expected 13.1`) unless `g++-11` is installed (`sudo apt install g++-11`).
2. **Boost C++ Headers** (`libboost-dev`, optional): Used if installed (`boost/algorithm/string.hpp`, `boost/container/stable_vector.hpp`); `incl/data_structs.hpp` provides standard C++17 fallbacks (`std::vector` and `boost::trim`/`trim_right`) via `__has_include` when Boost is absent.
3. **btllib (v1.4.9/1.4.10)**: Header files (`incl/btllib/`) and static library (`lib/libbtllib.a`) are vendored directly in the repository (**do not modify**).

### Common Commands
```bash
# Build optimized binary
make

# Build with debug symbols
make debug

# Run built-in sanity check on test/test.fa and test/test.fq
make run_test

# Remove compiled binary
make clean
```

---

## 6. Known Implementation Caveats (See `CODEREVIEW.md`)

Developers modifying or extending KoNuSeg should keep the following implementation details in mind (note: **all critical issues C1–C6 from `CODEREVIEW.md` have been resolved** without modifying `btllib`):
1. **Do Not Modify `btllib`**: Always respect `btllib` API contracts from the caller side (e.g., passing byte counts `(bfBits + 7) / 8` to `SeedBloomFilter` and `cbfSlots * sizeof(uint16_t)` to `CountingBloomFilter16`, and checking `SeedNtHash::roll()` / `!hit_seeds.empty()` on sequences that may contain `N`s).
2. **BED Boundary & Merging Quirks (`H1`, `H2`)**: `createSegments()` currently writes `segEndBase - 1` while advancing `segStartBase = segEndBase` (leaving 1-bp gaps between adjacent BED records under half-open BED semantics), and `mergeSegments()` still retains stale `segCpNo` on its final flush (`src/mergeSegments.cpp:124`).
3. **Global Mutable State**: `k_assembly`, `k_reads`, and `segTable` are defined as `inline` globals in `incl/data_structs.hpp`.
