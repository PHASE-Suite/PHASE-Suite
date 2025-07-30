# Analysis of an ecNGS Workflow Pipeline Bash Script

---

## Executive Summary
This report provides a comprehensive analysis of a bash script designed to execute the **PHASE** error-corrected Next-Generation Sequencing (ecNGS) workflow. The script systematically processes raw sequencing data (BAM files) through a series of bioinformatics steps, including alignment, quality filtering, identification of genomic regions, extraction of sequence mismatches, and aggregation of results. It leverages a suite of standard bioinformatics tools such as samtools, bwa mem, bedtools, jvarkit sam2tsv.jar, and GNU parallel to achieve its objectives. The pipeline is structured into modular functions, enabling parallel processing of multiple samples. The final output is a filtered and merged dataset of potential variants across all samples, optimized for downstream analysis in ecNGS applications.

---

## 1. Introduction to Error-Corrected Next-Generation Sequencing (ecNGS)
Next-generation sequencing (NGS) has revolutionized genomic analysis by enabling rapid and high-throughput sequencing of DNA and RNA. The general workflow for NGS involves several key stages: nucleic acid isolation, library preparation (fragmentation and adapter ligation), clonal amplification, sequencing by synthesis, and sophisticated bioinformatics data analysis. During data analysis, raw sequencing reads undergo base calling, alignment to a reference genome, and subsequent variant calling to identify differences between the sequenced DNA and the reference.

While standard NGS is powerful, its resolution for detecting low-frequency mutations can be limited by inherent sequencing errors and PCR artifacts. Error-corrected Next-Generation Sequencing (ecNGS) addresses this limitation by employing strategies, such as consensus sequencing, that compare reads derived from both original strands of DNA molecules. The increased sensitivity and accuracy of ecNGS make it a valuable tool for various applications, including detecting mutations in any cell type or tissue, quantifying drug- or chemical-induced mutations, characterizing mutational spectra, and identifying early markers of carcinogenesis. This script implements a pipeline designed to process data generated from such an ecNGS workflow, focusing on the computational analysis phase.

---

## 2. Overall Script Functionality
The provided bash script orchestrates a multi-step bioinformatics pipeline for processing Next-Generation Sequencing (NGS) data, specifically tailored for an error-corrected (ecNGS) workflow. The primary goal of the script is to take input BAM files (Binary Alignment Map format), align them to a reference genome, apply stringent quality and mapping filters, identify specific genomic regions, extract mismatch information, and then consolidate and further filter these results across all samples. The pipeline is designed for parallel execution of individual sample processing, enhancing computational efficiency.

The workflow can be broadly categorized into five main stages:
* **Initialization and Configuration**: Sets up the script's environment, defines global parameters, and creates necessary output directories.
* **Per-Sample Processing (Parallelized)**: For each input BAM file, it performs:
    * **Alignment and Initial Filtering**: Converts BAM to FASTQ, re-aligns reads, and filters them based on mapping quality and SAM flags.
    * **BED File Generation and Base Counting**: Converts filtered alignments to BED format, identifies overlapping regions from paired-end reads, and quantifies covered bases.
    * **Read Splitting and Mismatch Extraction**: Divides aligned reads into chunks and extracts mismatch information for both forward and reverse strands.
    * **Strand Joining and Sample-Specific Finalization**: Combines and sorts mismatch data from forward and reverse strands for each sample.
* **Global Result Merging and Final Filtering**: After all samples are processed, their individual results are concatenated and subjected to a final set of stringent filters to remove potential artifacts and duplicates, producing a refined dataset of high-confidence variants.
* **Error Handling**: A robust error trapping mechanism is in place to catch and report failures during script execution.
* **Timing and Reporting**: Tracks and reports the runtime for each sample and provides a summary of the overall process.

This modular design, coupled with parallel execution, is characteristic of modern bioinformatics pipelines, aiming for both efficiency and reproducibility in large-scale genomic data analysis.

---

## 3. Detailed Script Analysis

### 3.1. Initial Setup and Configuration
The script begins with a series of `set` commands and an error handler, which are crucial for robust shell scripting practices.

* `#!/usr/bin/env bash`: This shebang line ensures the script is executed with bash, regardless of the user's default shell.
* `set -euo pipefail`: This command activates several strict modes for the bash shell.
    * `set -e` (or `errexit`): Causes the script to exit immediately if any command fails (i.e., exits with a non-zero status). This prevents silent failures and propagation of errors, which is critical for complex pipelines where downstream steps depend on successful completion of upstream ones.
    * `set -u` (or `nounset`): Causes the script to exit immediately if an attempt is made to use an uninitialized variable. This helps catch typos and logical errors related to variable assignment early in development.
    * `set -o pipefail`: Ensures that the exit status of a pipeline is the exit status of the rightmost command that exited with a non-zero status. Without this, a pipeline would only report the exit status of the last command, potentially masking failures in earlier commands within the pipe.
* `IFS=$'\n\t'`: The Internal Field Separator (IFS) variable controls how bash performs word splitting. By default, IFS includes space, newline, and tab. Setting it to only newline and tab prevents word splitting on spaces, which is beneficial when iterating over lists of items (e.g., filenames) that might contain spaces, thereby avoiding unexpected parsing issues.
* `set -o errtrace` (or `set -E`): This option ensures that any `trap ERR` handler is inherited by shell functions, command substitutions, and commands executed in a subshell environment. Without `errtrace`, the `trap ERR` might not be triggered if a command fails inside a function or subshell, leading to unhandled errors.
* `error_handler()`: This function is defined to provide custom error messages. It captures the exit code (`$?`), line number (`BASH_LINENO`), and the function name (`FUNCNAME`) where the error occurred. This level of detail is invaluable for debugging and troubleshooting pipeline failures, significantly improving the script's maintainability and user-friendliness. The `trap error_handler ERR` command registers this function to be executed whenever a command exits with a non-zero status (i.e., an error occurs).

#### Configuration Variables
The script defines several key parameters at the beginning:

* `threads_per_job=8`: Specifies the number of CPU threads to be allocated per individual job (e.g., per sample processing).
* `jobs=3`: Dictates the maximum number of parallel `process_sample` jobs that GNU parallel will run concurrently.
* `splits=4`: Determines how many chunks each BAM file will be divided into for parallel processing within the `split_and_extract` function.
* `ref_fa="/path/to/reference/genome.fa"`: Path to the reference genome in FASTA format, essential for alignment and variant calling.
* `sam2tsv_jar="/path/to/jvarkit/dist/sam2tsv.jar"`: Path to the jvarkit `sam2tsv.jar` utility, used for converting SAM/BAM alignments to a tab-separated format for easier parsing.
* `workdir="$(pwd)"`: Sets the working directory to the current directory where the script is executed.
* `bam_dir="$workdir/original-bams"`: Specifies the directory containing the input BAM files.
* `out_dirs=(...)`: An array listing all necessary output subdirectories (bases, bed, output-bams, filtered-bams, fwd-rev, filtered-results, v1-joined). These directories are created using a `for` loop and `mkdir -p`, ensuring that all required output locations exist before processing begins.
* `summary_file="$workdir/sample_times.tsv"`: Path for a file to log the processing time for each sample. The `: > "$summary_file"` command truncates (or creates) this file at the start of the script.

### 3.2. Workflow Stages (Functions)
The core logic of the ecNGS pipeline is encapsulated within several bash functions, each performing a distinct step.

#### 3.2.1. Align & Filter (`align_and_filter`)
This function performs the initial alignment and filtering of reads for a given BAM file. It takes the input BAM filename as its first argument and extracts the sample name by removing the `.bam` extension.

The central command is a pipeline:
```bash
samtools fastq -@ 2 "$bam_dir/$bamfile" \
| bwa mem -t"$((threads_per_job-4))" -B21 -O11 "$ref_fa" - \
| tee >(samtools view -b -@ 1 -F 4079 -q 60 -o "$workdir/filtered-bams/${sample}.filtered.bam") \
| samtools view -b -@ 1 -o "$workdir/output-bams/${sample}.bam"
```

##### Command Breakdown
* **`samtools fastq`**: This command converts the input BAM file back into FASTQ format and streams it to standard output. This step effectively "un-aligns" the reads from their original BAM file, preparing them for re-alignment. This is often done in ecNGS workflows to ensure consistent processing from a raw-read-like state, especially if the original BAM might have undergone prior, less stringent processing.
* **`bwa mem`**: The FASTQ reads from `samtools fastq` are piped as standard input (`-`) to `bwa mem` for alignment against the reference genome (`$ref_fa`). `bwa mem` is a widely used aligner, recommended for high-quality reads typically 70-100bp or longer, and is known for its speed and accuracy.
    * `-t"$((threads_per_job-4))"`: Allocates a specific number of threads for `bwa mem`. For `threads_per_job=8`, this means `bwa mem` uses 4 threads. This thread allocation balances resource usage across the pipeline components.
    * `-B21`: This option sets the mismatch penalty to 21. A higher penalty makes the aligner less tolerant to mismatches, leading to more stringent alignment.
    * `-O11`: This option sets the gap open penalty to 11. A higher penalty discourages the introduction of gaps (insertions or deletions) in the alignment. These stringent penalties (`-B` and `-O`) are likely chosen to ensure only high-quality, well-matched reads are considered, which is crucial for ecNGS where distinguishing true mutations from sequencing errors is paramount.
* **`tee >(...) | samtools view`**: The output of `bwa mem` (in SAM format) is piped to `tee`. `tee` duplicates its standard input to multiple outputs.
    * The first output stream goes to a process substitution `>(...)`, which runs `samtools view` to create a **filtered** BAM file (`.filtered.bam`).
        * `-F 4079`: This is a critical filter based on SAM flags. The `-F` option excludes reads where any of the specified bits are set. 4079 in binary is `111111101111`. By using `-F 4079`, the script explicitly discards reads that are unmapped, have unmapped mates, are reverse complemented, are not primary alignments, fail quality checks, are duplicates, or are supplementary alignments. This is a very stringent filtering step, aiming to retain only high-confidence, uniquely mapped, non-duplicate, primary alignments on the forward strand, which is essential for accurate variant calling in ecNGS.
        * `-q 60`: Filters reads with a mapping quality (MAPQ) score less than 60. MAPQ scores quantify the probability that a read is misplaced, typically on a Phred scale, where a score of 60 indicates a very low probability of misalignment (1 in 1,000,000). This ensures only highly confident alignments are passed to the next stage, further reducing false positives in variant detection.
    * The second output stream from `tee` is piped to another `samtools view` command, which creates an **unfiltered** BAM file (`.bam`) in the `output-bams` directory. This file contains all aligned reads (before the stringent `-F 4079 -q 60` filters) and will be used later for mismatch extraction.

#### 3.2.2. BED-making & Base Counting (`make_bed_and_count`)
This function processes the `filtered.bam` file to generate a BED file and count the total number of bases covered.

```bash
bedtools bamtobed -i "$filt" \
| awk '...' OFS=$'\t' \
| tee "$workdir/bed/${sample}.bed" \
| awk '{ sum += ($4 - $3) } END { print sum }' > "$workdir/bases/${sample}.bases"
rm -f "$filt"
```

##### Command Breakdown
* **`bedtools bamtobed`**: Converts the filtered BAM alignment file into BED format. BED (Browser Extensible Data) format is a tab-delimited file used to define genomic regions.
* **`awk '...'`**: This `awk` script processes the BED output from `bedtools bamtobed`. Its purpose is to identify and extract overlapping regions from paired-end reads, specifically looking for "fwd" and "rev" read names. This logic identifies the common segment covered by both paired-end reads, which is crucial for ecNGS to ensure that mutations are supported by both strands, reducing false positives.
* **`tee "$workdir/bed/${sample}.bed"`**: The output of the `awk` script (the overlapping BED regions) is saved to a sample-specific BED file in the `bed` directory.
* **`awk '{ sum += ($4 - $3) } END { print sum }'`**: This `awk` command takes the output from the previous `tee` (the overlapping BED regions) and calculates the total length of these regions. It sums `(end - start)` for each line, effectively counting the total number of bases covered by the identified overlapping paired-end segments.
* **`rm -f "$filt"`**: The temporary `filtered.bam` file is removed to conserve disk space, as its processed information has now been captured in the BED and bases files.

#### 3.2.3. Split & Extract Mismatches (`split_and_extract`)
This function focuses on extracting mismatch information from the **unfiltered** aligned BAM file, parallelizing the process by splitting the BAM into chunks.

##### Command Breakdown
* **Header Extraction**: `samtools view -H "$bam" > "$hdr"` extracts the header from the BAM file and saves it to a temporary SAM header file.
* **Chunking**: `total=$(samtools view -c "$bam")` counts the total number of reads in the BAM file. `chunk=$(( (total + splits - 1) / splits ))` calculates the approximate number of reads per chunk. `samtools view "$bam" | split -l "$chunk" - "$prefix"` then splits the BAM file (without its header) into multiple smaller files.
* **Parallel Mismatch Extraction**: The script iterates over each generated chunk file. Inside the loop, for each chunk and for both `fwd` and `rev` strands, a subshell is launched in the background (`&`), enabling parallel execution.
    * `cat "$hdr" "$f"`: Concatenates the original SAM header with the current chunk file.
    * `java -jar "$sam2tsv_jar" -R "$ref_fa"`: The combined header and chunk data are piped to `jvarkit sam2tsv.jar`. This Java utility converts SAM alignments into a tab-delimited format, using the reference genome to determine reference bases at aligned positions.
    * `awk -v st="$st" '...'`: This `awk` script processes the `sam2tsv.jar` output to identify mismatches. It filters for lines where the read name ends with `/fwd` or `/rev` and the CIGAR operation is a mismatch. It checks if the read base is different from the reference base (`alt != ref`) and the reference is not 'N'. It then constructs a unique key for each mismatch: `readname:chromosome:reference_position:reference_base:alternative_base`.
* **`wait`**: After all background processes are launched, `wait` ensures that all parallel subshells complete before proceeding.
* **`rm -f "$hdr" ${prefix}*`**: Cleans up the temporary SAM header file and all the split chunk files to free up disk space.

#### 3.2.4. Join Strands & Finalize (`join_strands`)
This function aggregates and processes the mismatch data collected from forward and reverse strands for a single sample.

##### Command Breakdown
* **`cat` and `sort`**: The script concatenates all forward-strand mismatch files and all reverse-strand mismatch files into two separate `.all.fwd` and `.all.rev` files. It then sorts each file by the first field (the key), which is a prerequisite for the `join` command.
* **`join`**: `join -t $'\t' -j1` performs a database-style join operation on the first field of both files. This crucial command pairs up mismatch entries that share the exact same key, effectively identifying variants that are observed on **both** DNA strands of the original molecule. This is a core step in ecNGS, providing strong evidence against random, strand-specific sequencing errors.
* **Reformatting**: `tr ":" "\t"` replaces colons in the key field with tabs, splitting it into its constituent parts. An `awk` command then reorders the columns into a consistent output format.
* **Final Sort and Unique**: `sort -S80% --parallel=4 -u` sorts the reordered data and outputs only unique lines, removing duplicate entries.
* **Cleanup**: All intermediate files are removed to manage disk space.

### 3.3. Main Execution Flow
The script's main execution block orchestrates the per-sample processing and final aggregation.

* **`process_sample()`**: This wrapper function encapsulates the entire workflow for a single BAM file. It records the start time, calls the processing functions in sequence, and then calculates and prints the elapsed time for the sample.
* **`export`**: The script exports all defined functions and global variables, making them accessible to the subshells that GNU parallel executes.
* **`parallel --jobs "$jobs" process_sample ::: *.bam`**: This is the core parallelization command. It uses GNU Parallel to execute the `process_sample` function for each `.bam` file in the input directory, running up to the specified number of jobs concurrently.

### 3.4. 5) Combine & Final Filter
After all individual samples have been processed in parallel, this final stage merges their results and applies additional global filters.

* **`cat`**: All sample-specific `.FR.joined.v1` files are concatenated into a single `all.joined` file.
* **Final Filtering Pipeline**:
    * `sort -k6 -S80%`: Sorts the combined file by the 6th field (chromosome).
    * `uniq -u -f5`: This command filters for unique lines, ignoring the first 5 fields to determine uniqueness based on the variant key (chromosome, position, ref, alt). This step is a critical component of error correction, aiming to remove any remaining technical duplicates.
    * `awk '$2==0 || $2==16'` and `awk '$4==0 || $4==16'`: These filters ensure that the forward and reverse reads contributing to the variant were mapped in a standard orientation based on their SAM Flags.
    * `awk '$3==60 && $5==60'`: This filter enforces an extremely high mapping quality (MAPQ 60) for **both** the forward and reverse reads supporting the variant, dramatically reducing the likelihood of false positive calls.
* **Final Formatting and Deduplication**:
    * `cut` and `awk`: These commands select and reorder the columns to the standard format: `chromosome, position, reference_base, alternative_base, readname`.
    * `sort -k5 -S80% | uniq -u -f4`: This final deduplication step operates on the `readname` field (the 5th field after reordering). The `-u` option ensures that only lines with unique readnames are kept, effectively removing any remaining duplicate reads.

The final, highly filtered and formatted variant data is saved to `all.joined.2`.
