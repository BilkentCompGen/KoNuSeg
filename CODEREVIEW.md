# KoNuSeg Comprehensive Code & Manuscript Alignment Review

**Reviewed Repository**: `KoNuSeg` (`main` @ `b0e3231`)  
**Reference Manuscript**: Zambaku et al. (2026), *"Space-efficient detection of repeats in genome assemblies using Bloom filters"*  
**Scope**: All C++ sources (`src/*.cpp`), headers (`incl/*.hpp`), vendored `btllib` interface usage (`incl/btllib/*.hpp`), build system (`Makefile`), documentation (`README.md`), test fixtures (`test/*`), and consistency with the methods, equations, tables, and claims in the manuscript.

---

## 1. Executive Summary

KoNuSeg implements a clean, well-motivated algorithmic idea: using a two-tier spaced-seed **Bloom Filter (`SBF`) + Counting Bloom Filter (`CBF16`)** pipeline backed by `btllib` to filter sequencing errors in streaming reads and detect collapsed assembly duplications without read mapping. The codebase is compact (~580 lines of project C++ across 6 `.cpp` and 7 `.hpp` files) and straightforward to follow.

However, cross-checking the implementation line-by-line against `btllib` API contracts, C++17 memory-safety rules, and the attached manuscript revealed **6 Critical bugs**, **5 High-severity logic/IO issues**, **6 Manuscript-to-Code discrepancies**, and **major performance opportunities** (an estimated **5×–20× speedup** by using rolling `SeedNtHash` as intended rather than slicing `std::string::substr` per base). Several of these code behaviors also directly explain empirical anomalies discussed in the paper (such as elevated memory footprint on human data, residual hash-collision spikes, unique-sequence absorption during merging, and increased repeat calls at 10% read subsampling).

### Summary Severity Matrix

| ID | Severity | Location | Summary |
| :--- | :---: | :--- | :--- |
| **C1** | **Critical (Fixed)** | [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L86-L119) | Stack buffer overflow (`j <= windowSize`) in `noiseFilter()` + bug that disables outlier median filtering across almost every k-mer |
| **C2** | **Critical (Fixed)** | [loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L97-L103), [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L124-L138) | Out-of-bounds `hit_seeds[0]` vector read and uninitialized `repeats[0]` return when reads or contigs contain `N` / non-`ACGT` bases |
| **C3** | **Critical (Fixed)** | [populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L7-L26), [loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L12-L34) | Bit count $m$ passed to `btllib` constructors expecting **bytes**: **8×** over-allocation in `SBF` and **0.5×** counter slots (~15.8% FPR) in `CBF16` |
| **C4** | **Critical (Fixed)** | [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L21-L82) | Final segment of **every multi-segment contig** is silently dropped (`first_check == 0`), and the last k-mer of every contig is skipped |
| **C5** | **Critical (Fixed)** | [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L35-L45) | Division-by-zero crash (`SIGFPE`) when `k_reads / k_assembly < 0.65`, plus integer-division truncation of copy numbers |
| **C6** | **Critical (Fixed)** | [konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L5-L152) | `std::stoi` signed 32-bit overflow crashes (`std::out_of_range`) when `-anum` or `-knum` exceeds `2,147,483,647` (e.g., human genome size `3e9`) |
| **H1** | **High (Fixed)** | [mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L48-L143) | Stale `segCpNo` (set only on the 1st segment of the entire run) corrupts the final merged segment and contig-boundary checks; uninitialized vars on empty input |
| **H2** | **High (Fixed)** | [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L63-L95), [loadAssembly.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadAssembly.cpp#L28-L45) | 1-bp gaps (`segEndBase - 1`) between adjacent BED intervals leave unmasked bases in `bedtools maskfasta`; FASTA header comments corrupt BED chrom column |
| **H3** | **High (Fixed)** | [loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L4-L115), [konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L204-L280) | Custom FASTQ parser misparses Phred+33 quality lines starting with `@` (`Q=31`); missing file-open / `nullptr` checks (`SIGSEGV`/`SIGFPE`); no `.gz` support |
| **H4** | **High (Fixed)** | [konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L131-L141), [populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L26) | `SBF` initialized with `seed1.size()` instead of `kMerSize`; changing `-k` without `-s` silently produces an empty Bloom filter |
| **H5** | **High (Fixed)** | [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L63-L79), [mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L85-L89) | Tautological self-comparison `contigName.compare(contigName)` and unsigned integer underflows (`start - 2`, `start - segEndBase`) |
| **M1–M6** | **Medium (Fixed)** | Paper vs. Code | Quantized integer coverage at low depth, reversed Table 2 footnote / unused `-c_thresh` in Stage 2, dead `-g` flag, unexposed filter serialization |
| **P1–P3** | **Medium (Fixed)** | [populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L44-L95), [loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L11-L92), [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L126-L198) | Non-rolling `substr` + repeated heap allocations per k-mer (5×–20× slowdown), per-read OpenMP fork/join, data race on `s_pos`, non-deterministic output order |
| **L1–L5** | **Low** | Build / Hygiene | Memory leaks (`segTable` elements, `valSegment`), hardcoded `/usr/bin/g++-11`, undocumented/removable Boost dependency, tracked binaries, no `.gitignore` |

---

## 2. Critical Bugs (Memory Safety, Crashes, & Silent Data Corruption)

### C1. Stack Buffer Overflow & Disabled Median Outlier Filtering in `noiseFilter()`
- **Location**: [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L101-L133)
- **Paper Specification (p. 3 & Fig. 1d)**: *"To perform this filtering, we use a user-defined window size to slide over kmers and evaluate each k-mer based on the median count within that window."*
- **Code**:
  ```cpp
  uint slidingWindow[windowSize];          // Valid indices: 0 .. windowSize - 1
  for (uint i = 0; i < cSize; i++) {
      wIndex = 0;
      for (uint j = 0; j <= windowSize; j++) {   // BUG 1: runs windowSize + 1 times!
          if (i + j < cSize) {
              slidingWindow[j] = contigs->at(i + j).numRepeats; // BUG 2: OOB write at j == windowSize!
              wIndex++;
          }
      }
      if (wIndex == windowSize) {          // BUG 3: wIndex is windowSize + 1 for i < cSize - windowSize!
          sort(slidingWindow, slidingWindow + windowSize);
          uint median = slidingWindow[(windowSize / 2)];
          contigs->at(i + (windowSize / 2)).numRepeats = median; // BUG 4: mutates input in-place!
      }
  }
  ```
- **Impact**:
  1. **Undefined Behavior (Stack Buffer Overflow)**: For every k-mer `i < cSize - windowSize`, the inner loop executes `j == windowSize`, writing `4` bytes past the end of the stack VLA `slidingWindow`.
  2. **Outlier Filtering Practically Disabled**: For every k-mer `i < cSize - windowSize`, `wIndex` increments to `windowSize + 1` (e.g., `6` when `-w 5`). Consequently, `if (wIndex == windowSize)` (`6 == 5`) evaluates to **`false` across the entire contig**, executing **only once** at `i = cSize - windowSize`! Empirically, running `./bin/konuseg.out` on `test/test.fa` leaves isolated single-kmer hash-collision spikes (`seq3_1_cn: 436..468` with count `4` flanked by `2`s) completely unfiltered.
  3. **In-Place Recursive Mutation & VLA**: Even once `j < windowSize` is fixed, writing the median directly back into `contigs->at(i + windowSize / 2).numRepeats` causes subsequent windows (`i + 1`, ...) to compute medians over *already-filtered* values rather than raw CBF counts, and `uint slidingWindow[windowSize]` is a non-standard C++ VLA extension.
- **Fix**: Allocate a `std::vector<uint> window(windowSize)` and a separate `std::vector<uint> filtered(cSize)` (initialized from raw `numRepeats`), loop `for (uint j = 0; j < windowSize; ++j)`, and copy `filtered` back after scanning the contig.

---

### C2. Out-of-Bounds Vector Read & Uninitialized Return Value on Ambiguous Bases (`N`s)
- **Locations**: [src/loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L96-L98), [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L138-L152), [src/populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L44-L46)
- **Mechanism**: `btllib::NtHash` / `SeedNtHash` skips any k-mer containing non-`ACGT` characters (such as `N` gaps in draft contigs or ambiguous bases in reads; see [nthash.hpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/incl/btllib/nthash.hpp#L90-L97)). When passed a `k`-bp string containing `N`:
  1. **In `loadReads.cpp:96-98`**:
     ```cpp
     auto hit_seeds = seed_bf->contains(temp);
     if (find(hit_seeds[0].begin(), hit_seeds[0].end(), 0) != hit_seeds[0].end())
     ```
     `seed_bf->contains(temp)` returns an **empty** `std::vector<std::vector<unsigned>>` (`size() == 0`). Accessing `hit_seeds[0]` is an **out-of-bounds heap read (Undefined Behavior)** that either segfaults or reads stale heap pointers and falsely increments `k_reads`!
  2. **In `checkRepeats.cpp:138-152`**:
     ```cpp
     int repeats[seeds.size()]; // UNINITIALIZED!
     btllib::SeedNtHash h(seq, seeds, cbf->get_hash_num(), seeds[0].size());
     while (h.roll()) { ... repeats[i] = (unsigned)count; }
     return repeats[0];         // Returns stack garbage when h.roll() is false!
     ```
     When a contig contains `N`s (common in scaffolded/gap-filled drafts), `h.roll()` returns `false` immediately and `query()` returns **uninitialized stack memory** (typically the previous valid k-mer's count on that thread's stack), falsely annotating `N` gaps as repeats!
  3. **In `populateBF.cpp:44-46`**: `k_assembly` is unconditionally incremented for every `j` even when `temp` contains `N` and `seed_bf->insert()` inserts nothing, skewing `k_reads / k_assembly` if contigs contain `N`s.
- **Fix**: Check `!hit_seeds.empty()` before `hit_seeds[0]`; initialize `int repeats = 0;` in `query()`; and only count valid `ACGT` k-mers (where `SeedNtHash::roll()` returns `true`) toward `k_assembly`.

---

### C3. Bloom Filter & Counting Bloom Filter Sizing Unit Mismatch (Bits Passed as Bytes)
- **Locations**: [src/populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L10-L26), [src/loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L12-L32)
- **Mechanism**:
  The classic optimal Bloom filter size formula computes $m$ in **bits** (or **slots/counters** for a Counting Bloom filter):
  $$m = -\frac{n \ln p}{(\ln 2)^2}$$
  However, both [btllib::SeedBloomFilter](file:///home/calkan/code/bilkentcompgen/KoNuSeg/incl/btllib/bloom_filter.hpp#L420-L428) and [btllib::CountingBloomFilter16](file:///home/calkan/code/bilkentcompgen/KoNuSeg/incl/btllib/counting_bloom_filter.hpp#L48-L54) take **`size_t bytes`** (memory budget in bytes) as their first constructor argument!
  1. **Seed Bloom Filter (`populateBF.cpp:26`) — 8× Memory Over-Allocation**:
     - `SBF *seed_bf = new SBF(bfSize, ...)` passes `bfSize` (bits) as `bytes`, allocating **`8 * bfSize` bits**!
     - Meanwhile, line 22 prints `Standard BF Memory: bfSize / (8 * 1024^3) GB` — **underreporting actual SBF RAM usage by 8×**!
     - *Empirical verification on `make run_test`*: `a_num = 1174, fpr = 0.01` gives `bfSize = 11252` bits. Optimal occupancy at 1,174 distinct k-mers should be ~`0.50`, yet stdout reports `BF occupancy: 0.0910959` ($1 - e^{-7 \times 1174 / (11252 \times 8)} \approx 0.087$). On human defaults (`-anum 1.5e9`), `SBF` prints `1.67 GB` while actually allocating **13.39 GiB (14.38 GB)**!
  2. **Counting Bloom Filter (`loadReads.cpp:32`) — 0.5× Counter Slots & 15.8× Higher FPR**:
     - `CBF` is aliased to `btllib::CountingBloomFilter16` (`uint16_t`, **2 bytes per counter**, [data_structs.hpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/incl/data_structs.hpp#L241)), contradicting the comment on line 24 (`// CBF uses 1 byte per slot`) and line 22 (`CBF Size: ... bits`).
     - Inside `CountingBloomFilter<uint16_t>`, `array_size = bytes / sizeof(uint16_t) = cbfSize / 2` ([counting_bloom_filter.hpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/incl/btllib/counting_bloom_filter.hpp#L663)).
     - Passing `cbfSize` as `bytes` gives the filter **only half the required counter slots (`m / 2`)**, raising the effective false-positive / hash-collision rate at full capacity from the target **1.0%** (`0.01`) to $(1 - e^{-2 \ln 2})^7 \approx$ **15.8%**!
- **Fix**:
  - For `SBF`: pass `size_t bfBytes = (bfBits + CHAR_BIT - 1) / CHAR_BIT;`.
  - For `CBF` (`CountingBloomFilter16`): pass `size_t cbfBytes = cbfSlots * sizeof(uint16_t);` and update the memory printout to divide `cbfBytes` by `1024.0^3` (with label `counters` instead of `bits`). Note: Fixing `SBF` saves ~11.7 GiB on human while fixing `CBF16` doubles its counter array from ~13.4 GiB to ~26.8 GiB (net +15 GiB, or switch to `CountingBloomFilter8` if counts `<= 255` suffice, saving 11.7 GiB net while restoring true 1% FPR!).

---

### C4. Final Segment Dropped on Every Multi-Segment Contig & Final K-mer Skipped
- **Location**: [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L45-L95)
- **Code**:
  ```cpp
  for (uint i = 1; i <= Contigs->size(); i++) {
      if (i < Contigs->size() - 1) {
          if (copyNumberInsideThreshold(...)) {
              ...
          } else {
              first_check = 0; // Set to 0 permanently after the first segment break!
              ... // Emits previous segment, resets segStartBase..segEndBase for new segment
          }
      }
  }
  if (first_check == 1) { // FALSE whenever >= 1 segment break occurred in the contig!
      valSegment->push_back(...);
      *wout << ...;
  }
  ```
- **Impact**:
  1. **Trailing Segment Lost on Every Multi-Segment Contig**: As soon as a contig experiences a single copy-number transition (`else` branch), `first_check` becomes `0`. When the loop finishes, the final accumulated segment extending to the end of the contig is **never flushed** to `valSegment` or `wout`!
     - *Empirical verification*: In `test/test.fa`, `seq3_1_cn` is **796 bp** long, yet in `test_result.bed` its output terminates at **788 bp** — the final segment at the 3' end of the contig was silently lost. On real chromosomes or contigs, the final unique or repetitive segment (potentially tens or hundreds of kilobases!) at the end of every multi-segment contig is completely missing from both unmerged and merged BED outputs.
  2. **Final K-mer Skipped**: `if (i < Contigs->size() - 1)` skips index `Contigs->size() - 1` (the last k-mer of every contig). Even single-segment contigs (`seq1_1_cn`: 304 bp, `seq2_1_cn`: 220 bp) stop 1 k-mer early (`302` and `218` after subtracting an additional `1` via `segEndBase - 1`).
- **Fix**: Iterate `for (size_t i = 1; i < Contigs->size(); ++i)` and **unconditionally** flush the active trailing segment (`segStartBase` to `segEndBase`) after the loop completes (removing `first_check`).

---

### C5. Division-by-Zero Crash (`SIGFPE`) & Integer Division Truncation in `createSegments()`
- **Location**: [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L33-L58)
- **Code**:
  ```cpp
  uint coverage_depth = roundCount(static_cast<float>(k_reads) / k_assembly);
  ...
  wMedian = Contigs->at(0).numRepeats / coverage_depth;
  ...
  if (copyNumberInsideThreshold(Contigs->at(i).numRepeats / coverage_depth, wMedian/kCount, count_threshold) ...)
  ```
- **Impact**:
  1. **`SIGFPE` Crash**: If `k_assembly == 0` (e.g., all contigs `< k`) OR if `float(k_reads) / k_assembly < 0.65` (ultra-low sequencing depth, high read error rate, mismatched FASTA/FASTQ, or empty/unreadable FASTQ), `roundCount()` returns `0.0f` or `0.5f`, which truncates to **`uint coverage_depth = 0`**, crashing immediately with **hardware integer division by zero (`SIGFPE`)** at line 40!
  2. **Integer Division Truncation**: Because both `Contigs->at(i).numRepeats` and `coverage_depth` have type `uint`, `Contigs->at(i).numRepeats / coverage_depth` on lines 40 and 50 performs **integer division** (e.g., `19 / 10 == 1` instead of `1.9f`) before converting to `float`! Meanwhile, lines 57 and 84 compute `(Contigs->at(i).numRepeats * 1.0 / coverage_depth * 1.0)` in floating point, creating an inconsistency between threshold checking and segment mean accumulation. (Note also that `wMedian` is a running *sum*, so `wMedian / kCount` is the running *mean*, not median.)
- **Fix**: Keep `float coverage_depth = static_cast<float>(k_reads) / static_cast<float>(k_assembly);` (guarded against `k_assembly == 0` and `coverage_depth <= 0.0f`, e.g., clamping or aborting cleanly with an informative error if `k_reads == 0`), and perform floating-point division `static_cast<float>(numRepeats) / coverage_depth` everywhere.

---

### C6. Signed 32-Bit Integer Overflow (`std::stoi`) on `-anum` and `-knum`
- **Location**: [src/konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L22-L23)
- **Code**:
  ```cpp
  else if (arg == "-knum" && i + 1 < argc) config->k_num = stoi(argv[++i]);
  else if (arg == "-anum" && i + 1 < argc) config->a_num = stoi(argv[++i]);
  ```
- **Impact**: `std::stoi` parses a signed 32-bit `int` (max `2,147,483,647`). Both `k_num` and `a_num` are `ulong`, and [README.md](file:///home/calkan/code/bilkentcompgen/KoNuSeg/README.md#L75) advises: *"If `-anum` and `-knum` are not known, the size of the genome can be used"* (where human genome size `-g` is `3,000,000,000`). Passing `-anum 2500000000` or `-knum 3000000000` immediately crashes KoNuSeg with an uncaught `std::out_of_range` exception (`SIGABRT`). Furthermore, negative values passed to `stoi` for `uint` parameters (`-k`, `-w`, `-minseg`, `-c_thresh`) silently wrap around to `~4.29e9`.
- **Fix**: Use `std::stoull` for `ulong`/`uint64_t` counts, validate parameter ranges (`k > 0`, `0 < fpr < 1`, `w >= 1`, `seed1.size() == k_len`), and wrap CLI numeric parsing in `try / catch (const std::exception&)`.

---

## 3. High-Severity Logic & Input/Output Bugs

### H1. Stale `segCpNo` & Uninitialized Variables in `mergeSegments()`
- **Location**: [src/mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L53-L126)
- **Issues**:
  1. `segCpNo` is assigned **only once** at line 78 (`if (currentSegment == "")`) on the **very first segment of the entire run** and is never updated again!
     - At line 110 (contig transition): `if (roundCount(segCpNo) > 0)` checks the *first contig's first segment*, not the segment currently being flushed! If the first segment of the run had rounded copy number `0`, **the final merged segment of every contig in the assembly is silently dropped**.
     - At line 125 (final segment after loop): `valSegment->push_back(SegmentInfo(currentSegment, segStartBase, segEndBase, roundCount(segCpNo)));` ignores `processData(values)` completely and overwrites the last segment's copy number with that of the very first segment in `segTable`! Because OpenMP populates `segTable` in non-deterministic thread completion order, this value also varies randomly across runs.
  2. `segStartBase`, `segEndBase`, and `segCpNo` (lines 53–55) are uninitialized; if `segTable` has no segments, line 125 writes uninitialized stack values to `<outfile>_merged.bed`.
- **Fix**: Guard `if (!values.empty())` after the loop and flush using `processData(values)`; remove `segCpNo` completely (or check `processData(values) > 0` consistently).

### H2. 1-bp BED Coordinate Gaps & FASTA Header Descriptions Corrupting BED Column 1
- **Locations**: [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L65-L77), [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L192), [src/loadAssembly.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadAssembly.cpp#L31-L36)
- **Issues**:
  1. **1-bp Gaps Between Adjacent Segments**: BED format (used by `bedtools intersect` and `bedtools maskfasta` in the paper's evaluation, pp. 3 & 11) uses **0-based, half-open `[start, end)`** coordinates. In `createSegments()`, k-mer `j` has `[j, j + k_len)`. Outputting `segEndBase - 1` as column 3 while setting the next segment's `segStartBase = segEndBase` produces intervals like `0\t32`, `33\t66`, `67\t105` — **skipping 1 base (`32`, `66`, ...) at every single segment boundary**! When passed to `bedtools maskfasta`, those boundary bases remain unmasked (`A/C/G/T` amidst `N`s).
  2. **No-Op `line.erase(0, 0)` & Unsplit FASTA Headers**:
     In [loadAssembly.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadAssembly.cpp#L32), `line.erase(0, 0);` erases **0 characters** (`erase(pos, count)`), leaving `'>'` attached to `Sequence::name` (worked around later via `name.substr(1, ...)` in `checkRepeats.cpp:192`). Worse, if a FASTA header contains standard metadata after whitespace (e.g., `>ctg1 len=54200`), the entire string `ctg1 len=54200` is emitted into BED column 1, breaking tab/whitespace-delimited BED compatibility with `bedtools`!
- **Fix**: Strip `>` properly (`line.substr(1)`), truncate `Sequence::name` at the first whitespace (`line.substr(0, line.find_first_of(" \t"))`), and emit true half-open `[segStartBase, segEndBase)` coordinates without `- 1`.

### H3. Fragile Custom FASTQ/FASTA Parsers & Missing File / Null Checks
- **Locations**: [src/loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L34-L61), [src/konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L85), [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L174), [src/mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L51)
- **Issues**:
  1. **Phred+33 `@` Quality Score Bug**: In `loadReads.cpp:41`, `if (line[0] == '@')` is used to detect FASTQ headers without tracking line state (4-line FASTQ structure or sequence vs. quality length). In Phred+33 encoding, **`'@'` (ASCII 64) is a valid, common quality score (`Q = 31`)**! On multi-line FASTQ files (which lines 56–60 attempt to support), a quality line starting with `@` sets `read_counter = 1`, causing subsequent wrapped quality lines to be concatenated into `content` and **hashed into the Counting Bloom Filter as DNA sequence**!
  2. **Unchecked File Streams & Null Dereference**:
     - `konuseg.cpp:85`: If `loadAssembly()` fails to open `-fa` and returns `nullptr`, `main()` does not check for `nullptr` and immediately crashes with **`SIGSEGV`** at `populateBF.cpp:33` (`assembly->size()`).
     - `loadReads.cpp:34`: `ifstream reads(readsFile);` never checks `reads.is_open()`; a typo in `-fq` silently reads 0 k-mers and crashes with **`SIGFPE`** in `createSegments()`.
     - `checkRepeats.cpp:174` & `mergeSegments.cpp:51`: `ofstream wout(...)` is never checked for write errors.
  3. **No Compressed Input Support (`.gz`)**: Sequencing reads (`374 GB` in human!) are almost universally stored as `.fastq.gz` or `.fasta.gz` (and FASTA reads are also common, e.g., `long_reads.fa` in Supplementary Protocol p. 10–11). `btllib::SeqReader` is **already compiled into `lib/libbtllib.a`** and handles multi-threaded parsing of FASTA, FASTQ, `.gz`, `.bz2`, `.xz`, and SAM/BAM automatically!
- **Fix**: Replace both custom `std::getline` parsers (`loadAssembly.cpp` and `loadReads.cpp`) with `btllib::SeqReader`, and check return status / output streams before proceeding.

### H4. Silent Failure When User Overrides `-k` Without `-s`
- **Location**: [src/populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L26)
- **Code**: `SBF *seed_bf = new SBF(bfSize, seed1.size(), {seed1}, optHashFunctions);`
- **Impact**: If a user runs `-k 21` without supplying a 21-character `-s` seed, `k_len` becomes `21` while `seed1.size()` stays `32`. `populateBF()` then extracts 21-bp substrings (`substr(j, 21)`) and passes them to a `SeedBloomFilter` expecting `k = 32`. `SeedNtHash::init()` immediately returns `false` (`k > seq_len`) for every k-mer, resulting in `BF occupancy: 0` and empty output with **zero warning or error**! Similarly, if `-s` is shorter than `-k`, `populateBF` hashes multiple rolling windows per `substr(j, kMerSize)`.
- **Fix**: Validate at startup in `parseCommandLine()`: if `-k` is provided without `-s`, either default `-s` to `std::string(k_len, '1')` (exact k-mers) or require a matching `-s`, and error out if `config->seed1.size() != config->k_len` or contains characters other than `'0'` and `'1'`.

### H5. Tautological String Comparison & Unsigned Integer Underflows
- **Locations**: [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L63-L65), [src/mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L89)
- **Issues**:
  1. `checkRepeats.cpp:63`: `Contigs->at(i).contigName.compare(Contigs->at(i).contigName) == 0` compares a string **to itself** (always `true`) instead of comparing against `currentContig`.
  2. `checkRepeats.cpp:65`: `(Contigs->at(i+1).start - 2) > segStartBase` — both operands are `uint`; if `start < 2`, `start - 2` wraps around to `4294967294`. Also, because of the `i < Contigs->size() - 1` check, `i + 1` stays in-bounds today, but fixing C4 (`i < Contigs->size()`) would make `Contigs->at(i + 1)` throw `std::out_of_range` unless rewritten cleanly!
  3. `mergeSegments.cpp:89`: `(segItem->at(i).start - segEndBase) <= 3 * k_len` — unsigned subtraction wraps around to `~UINT_MAX` if `segItem->at(i).start < segEndBase`.

---

## 4. Manuscript-to-Code Discrepancies & Scientific Insights

| # | Manuscript Section / Claim | Actual Code Behavior | Scientific / Practical Implication |
| :- | :--- | :--- | :--- |
| **M1** | **Methods (pp. 2–3) & Results (pp. 4–5, Fig. S3–S4)**: Each k-mer's count is divided by Bloom Filter Coverage (`k_reads / k_assembly`). At 10% subsampling in *S. lycopersicum* and *H. sapiens*, estimated repeat percentage unexpectedly *increased*. | `loadReads.cpp:75` prints `std::round(float(k_reads)/k_assembly)`, whereas `checkRepeats.cpp:33` uses `uint coverage_depth = roundCount(float(k_reads)/k_assembly)` (quantized to integer!). | At 10% ONT depth, effective CBF coverage is ~`1–2×` (Fig. S3–S4). Because `roundCount(1.64f)` rounds **down** to `1`, dividing counts by `1` instead of `1.64` inflates all copy estimates by **+64%**, pushing neutral `1.0×` regions above `1.65` (`roundCount -> 2`, duplicated)! Using exact `float` coverage eliminates this low-depth quantization bias. |
| **M2** | **Methods (p. 3), Table 2 Footnote (p. 4), & Table S1 (p. 9)**: Table S1 & CLI help define `-c_thresh` as *"Count threshold for merging"*. Table 2's footnote states *Unmerged* uses minimum segment size without merging similar copy estimates, while *Merged* merges segments *"according to their copy estimate values"*. | `-c_thresh` is **unused** in `mergeSegments()` (`src/mergeSegments.cpp:49`) and is actually the Stage 1 k-mer extension threshold in `createSegments()`. Conversely, `min_seg_len` (`-minseg`) is passed to `createSegments()` where it is **unused**, and is only used in `mergeSegments()`. Furthermore, `mergeSegments()` concatenates adjacent segments `< min_seg_len` **regardless of copy-number difference**! | Table 2's footnote reverses the roles of `minseg` and copy-estimate similarity. In `mergeSegments()`, pooling a 200 bp single-copy (`cn=1`) segment with a 400 bp repeat (`cn=10`) merges them unconditionally; and because `remove_count = size / 4` (`mergeSegments.cpp:24`) equals **`0`** whenever `< 4` segments are pooled, the high value is not discarded (`(1+10)/2 -> 5`), directly explaining the *"unique sequence absorption during segment merging"* noted in Discussion (p. 7). |
| **M3** | **Table S1 (p. 9) & `README.md`**: Lists `-g` (`Genome size, default: 3,000,000,000`) as an optional parameter. | Parsed into `Config::genome_size` (`konuseg.cpp:15`) and passed to `checkRepeats(..., long double genomeSize, ...)` (`checkRepeats.cpp:161`), where it is **completely unused**. | Users setting `-g 100000000` (without setting `-anum` and `-knum`) may assume KoNuSeg sizes its Bloom filters for a 100 Mbp genome, when it actually still allocates filters for the `1.5e9` default! If `-anum`/`-knum` are omitted while `-g` is supplied, `-g` could sensibly serve as the fallback for `a_num`/`k_num` (as suggested in `README.md:75`). |
| **M4** | **Methods (p. 2) & Fig. 2 (p. 5)**: Attributes slight peak-memory increase at higher thread counts to *"a larger volume of sequence reads to be retained in memory and processed simultaneously"*. | `loadReads.cpp:38-61` reads one FASTQ record at a time sequentially on the main thread (`string content`) before entering `#pragma omp parallel for` per read. | Only **one read** is ever retained in RAM regardless of thread count; the ~1 GB memory growth from 1 to 100 threads in Fig. 2 comes from glibc per-thread malloc arenas, OpenMP thread stacks, and concurrent `CONTIGLIST` / `std::string` allocations in `checkRepeats()`. |
| **M5** | **Discussion (p. 7)**: *"the ability to serialize filter states makes the tool well-suited for iterative research workflows, in which prebuilt filter index files can be reused"*. | `btllib` supports `.save(path)` and file-path constructors on both `SeedBloomFilter` and `CountingBloomFilter16`, **but KoNuSeg exposes no CLI options** to save or load `.bf` / `.cbf` files. | Adding `--save-cbf <path>` / `--load-cbf <path>` is ~15 lines of code using `btllib::CountingBloomFilter16::is_bloom_file()` and `.save()`, immediately fulfilling the manuscript's claim. |
| **M6** | **High-Coverage Copy Estimation (Table S2)**: *C. elegans* ONT dataset has **683×** coverage. | `CBF` uses 16-bit counters (`CountingBloomFilter16`, max `65,535`). | `btllib` safely saturates at `65,535` without wrap-around, which caps measurable copy number at $65,535 / \text{CBF\_cov}$ (~120 copies at 545× CBF coverage in *C. elegans*). Benign for binary `>1` / `>10` masking, worth documenting in `README.md`. |

---

## 5. Performance & Concurrency Bottlenecks

### P1. Non-Rolling `substr()` Slicing & Heap Allocations Inside Inner Loops (5×–20× Speedup Opportunity)
- **Locations**: [src/populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L42-L47), [src/loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L92-L119), [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L138-L194)
- **Issue**: `btllib::SeedNtHash` computes rolling spaced-seed hashes in **$O(1)$ time per base** across an entire sequence. However, all three stages slice `seq.substr(j, kMerSize)` at every base position `j`, forcing `SeedNtHash` to re-initialize from scratch in $O(k \cdot h)$ time on a 1-kmer string and never roll! Worse, for **every single k-mer** in `loadReads.cpp` (e.g., ~200 billion k-mers in the 68× human ONT dataset!):
  1. `temp = read.substr(j, kMerSize)` allocates/copies `k` bytes.
  2. `seed_bf->contains(temp)` constructs a `SeedNtHash` (heap-allocating internal hash/block arrays) and returns a `std::vector<std::vector<unsigned>>` (2 heap allocations).
  3. On a hit, `populate(temp, seeds, cbf)` constructs **a second identical `SeedNtHash`** on `temp`, heap-allocates `std::unique_ptr<uint64_t[]>(new uint64_t[hash_num])`, `std::copy`s the hashes (`cbf->insert(h.hashes())` accepts `const uint64_t*` directly without copying!), and inserts into `cbf`!
- **Recommendation**:
  - In `populateBF()`: simply call `seed_bf->insert(contig.data)` once per contig (or per parallel chunk of a contig with `k - 1` overlap)!
  - In `loadReads()` and `checkRepeats()`: construct **one** `btllib::SeedNtHash h(seq, parsed_seeds, hash_num, kMerSize)` per read/contig (pre-parsing `parsed_seeds` once via `btllib::parse_seeds`), loop `while (h.roll())`, check `if (seed_bf->contains(h.hashes()))`, and immediately call `cbf->insert(h.hashes())` — **reusing the exact same `h.hashes()` pointer with zero heap allocations and zero redundant hashing!** This alone can cut the 6-hour human runtime down dramatically (likely under 30–60 minutes).

### P2. Single-Threaded File I/O & Per-Read OpenMP Fork/Join
- **Location**: [src/loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L38-L61) & [L91](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L91)
- **Issue**: Reading a 374 GB FASTQ via unbuffered `std::getline` on 1 thread serializes all disk I/O, while spawning `#pragma omp parallel for` inside `populateCBF()` incurs an OpenMP barrier synchronization per read (severe overhead when short reads are used, which p. 2 notes are supported).
- **Recommendation**: Use `btllib::SeqReader` (`SeqReader::Flag::LONG_MODE` or `SHORT_MODE`) inside a single outer `#pragma omp parallel` block so each OpenMP worker pulls batches of reads concurrently and rolls `SeedNtHash` locally with zero synchronization (`btllib` filters use lock-free `std::atomic` operations internally).

### P3. OpenMP Data Race on `s_pos` & Non-Deterministic Output Order
- **Locations**: [src/populateBF.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L8) & [L45](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/populateBF.cpp#L45), [src/checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L24) & [L178](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L178)
- **Issues**:
  1. **Data Race (UB)**: In `populateBF.cpp:8`, `uint s_pos = 0;` is declared outside `#pragma omp parallel for` (making it implicitly `shared`) and written concurrently by all threads at line 45 (`s_pos = j + kMerSize;`) without synchronization. `s_pos` is never read and should be deleted.
  2. **Schedule Imbalance & Non-Deterministic BED Order**: `checkRepeats.cpp:178` uses default static loop scheduling over contigs (which can suffer severe load imbalance when chromosome-scale scaffolds mix with short contigs; `schedule(dynamic, 1)` is much better) and writes directly to `wout` / `segTable` in thread-completion order inside `#pragma omp critical`. Pre-sizing `segTable->resize(num_elements)` so thread processing contig `idx` writes to `(*segTable)[idx]` removes the critical section entirely and guarantees **100% deterministic, FASTA-ordered BED output**.

---

## 6. Code Quality, Memory Management, & Build Hygiene

- **L1. Memory Leaks**:
  - [src/konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L133): `delete segTable;` frees the outer vector but leaks all per-contig `SEGMENTINFO*` heap allocations created at [checkRepeats.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/checkRepeats.cpp#L42).
  - [src/mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L63): `SEGMENTINFO *valSegment = new SEGMENTINFO();` is never `delete`d. Using value-typed `std::vector<SegmentInfo>` eliminates manual `new`/`delete` across the entire project.
- **L2. Header Pollution & Global Mutable State ([incl/data_structs.hpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/incl/data_structs.hpp))**:
  - `using namespace std;` at global scope in the root header pollutes all translation units.
  - Mutable globals (`inline ulong k_assembly = 0;`, `inline ulong k_reads = 0;`, `inline SEGMENTTABLE* segTable = new SEGMENTTABLE();`) create hidden coupling across pipeline stages.
  - Missing standard headers (`<iomanip>` for `std::setprecision` in `loadReads.cpp:25`, `<algorithm>` for `std::sort`/`std::find`, `<fstream>`) rely on transitive includes from Boost/`btllib`.
- **L3. Build System, GCC 11 Slim-LTO Lock-In, & Boost Dependency ([Makefile](file:///home/calkan/code/bilkentcompgen/KoNuSeg/Makefile), [README.md](file:///home/calkan/code/bilkentcompgen/KoNuSeg/README.md))**:
  - `README.md` states *"GCC / G++ version **11** or higher"*, whereas `lib/libbtllib.a` was compiled with **GCC 11 slim LTO** (LTO bytecode version 11.3, zero-size `.text` sections). Consequently, GCC 12/13+ cannot link `lib/libbtllib.a` (`lto1: bytecode stream ... generated with LTO version 11.3 instead of the expected 13.1`); **exact `g++-11` is required** unless `libbtllib.a` is built without slim LTO (or with `-ffat-lto-objects`).
  - `incl/data_structs.hpp` previously required six Boost headers (four completely unused) that were undocumented in `README.md`; C++17 `__has_include` fallbacks have now been added in `incl/data_structs.hpp` so Boost is optional.
- **L4. Repository Hygiene**:
  - Compiled binary `bin/konuseg.out`, 8 MB static library `lib/libbtllib.a`, duplicate `cpptoml.h` copies, and `btllib`'s upstream test suite (`incl/tests/*`) are committed to Git, and there is no `.gitignore`.
- **L5. Minor CLI / Logging Polish**:
  - [src/konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L93) logs `Loading reads from: <fastq>` right before `Populating Seed Bloom Filter...` (copy-paste from Stage 2; Stage 1 processes assembly contigs).
  - [src/loadReads.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/loadReads.cpp#L75): typo `BLOOM FITLER COVERAGE DEPTH:`.
  - [src/mergeSegments.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/mergeSegments.cpp#L51): `outfile_path + "_merged.bed"` produces `foo.bed_merged.bed` when `-o foo.bed` is given (stripping a trailing `.bed` before appending `_merged.bed` is cleaner).
  - [src/konuseg.cpp](file:///home/calkan/code/bilkentcompgen/KoNuSeg/src/konuseg.cpp#L67): `return -1;` from `main()` exits with status `255` on POSIX; prefer `return 1;` (`EXIT_FAILURE`) and support `-h` / `--help` with exit code `0`.

---

## 7. Recommended Remediation Priority

1. **P0 — Fix Memory Safety, Sizing, & Data Loss Bugs (`C1`–`C6`, `H1`–`H2`)**:
   - Fix `noiseFilter()` bounds (`j < windowSize`), non-VLA storage, and out-of-place median output (**C1**).
   - Guard `N`-containing k-mers in `loadReads()`, `query()`, and `populateBF()` (**C2**).
   - Convert `bfSize` bits to bytes (`(bfSize + 7) / 8`) for `SBF` and `cbfSize * sizeof(uint16_t)` bytes for `CBF16` (**C3**).
   - Fix `createSegments()` loop bounds (`i < Contigs->size()`), always flush the final contig segment, and output standard half-open `[start, end)` BED coordinates (**C4**, **H2**).
   - Keep `coverage_depth` as `float` (guarding `k_assembly == 0` / `k_reads == 0`) (**C5**, **M1**).
   - Switch CLI parsing to `std::stoull` with validation (`seed1.size() == k_len`, etc.) (**C6**, **H4**) and fix `mergeSegments()` final-segment `processData(values)` call (**H1**).
2. **P1 — Replace Custom Parsers & Align CLI/Paper Semantics (`H3`, `M2`–`M5`)**:
   - Adopt `btllib::SeqReader` for multi-threaded FASTA/FASTQ (`.gz`/`.xz`) reading and clean header ID extraction (**H2**, **H3**).
   - Clarify or wire `-g` (e.g., default `a_num`/`k_num` to `genome_size` when unspecified) and `-c_thresh` (optionally enforce copy-state compatibility during Stage 2 merging to prevent unique-sequence absorption) (**M2**, **M3**).
   - Expose optional `--save-cbf` / `--load-cbf` flags via `btllib` serialization (**M5**).
3. **P2 — Unlock Rolling-Hash Performance & Clean Up Build (`P1`–`P3`, `L1`–`L5`)**:
   - Replace per-position `substr()` + `SeedNtHash` reconstruction with a single rolling `SeedNtHash` per sequence and pre-indexed deterministic contig output (**P1**–**P3**).
   - Remove unused Boost headers, generalize `Makefile` (`CXX ?= g++ -Wall -Wextra`), and add a `.gitignore` (**L1**–**L5**).
