# Technical Analysis of the PHASE Suite Workflow Pipeline

---

## Executive Summary

This report provides a comprehensive analysis of the **PHASE Suite**, a modular bioinformatics workflow designed for error-corrected Next-Generation Sequencing (ecNGS). The suite is composed of a primary Bash script (`phase_pipeline.sh`) and external configuration files that systematically process raw sequencing data (BAM files) through a series of bioinformatics steps. The pipeline performs alignment, quality filtering, consensus region identification, sequence mismatch extraction, result aggregation, and barcode-based demultiplexing. It leverages a suite of standard bioinformatics tools such as `samtools`, `bwa`, `bedtools`, `seqkit`, `sam2tsv`, and `GNU Parallel`. The final output is a set of demultiplexed, filtered data files optimized for direct input into the suite's interactive visualization application.

---

## 1. Introduction to Error-Corrected Next-Generation Sequencing (ecNGS)

Next-generation sequencing (NGS) has revolutionized genomic analysis by enabling rapid and high-throughput sequencing of DNA and RNA. The general workflow for NGS involves several key stages: nucleic acid isolation, library preparation (fragmentation and adapter ligation), clonal amplification, sequencing by synthesis, and sophisticated bioinformatics data analysis. During data analysis, raw sequencing reads undergo base calling, alignment to a reference genome, and subsequent variant calling to identify differences between the sequenced DNA and the reference.

While standard NGS is powerful, its resolution for detecting low-frequency mutations can be limited by inherent sequencing errors and PCR artifacts. Error-corrected Next-Generation Sequencing (ecNGS) addresses this limitation by employing strategies, such as consensus sequencing, that compare reads derived from both original strands of DNA molecules. The increased sensitivity and accuracy of ecNGS make it a valuable tool for various applications, including detecting mutations in any cell type or tissue, quantifying drug- or chemical-induced mutations, characterizing mutational spectra, and identifying early markers of carcinogenesis. This script implements a pipeline designed to process data generated from such an ecNGS workflow, focusing on the computational analysis phase.

---

## 2. Overall Script Functionality

The **PHASE Suite** orchestrates a multi-stage bioinformatics pipeline designed for a complete ecNGS workflow, from raw data to demultiplexed results. The suite's architecture emphasizes modularity, safety, and efficiency by separating the core logic (`phase_pipeline.sh`) from user-defined settings (`pipeline.conf`, `my_barcodes.cfg`).

The workflow is executed in three primary stages:

* **Part 1: Per-Sample Core Mutation Analysis (Parallelized):** Each input BAM file is processed independently through a four-step pipeline to produce a high-confidence, error-corrected list of mutations for that sample.
* **Part 2: Global Result Merging and Final Filtering:** After all samples are processed, their individual results are concatenated and subjected to a final set of stringent global filters to produce a master list of high-confidence variants.
* **Part 3: Per-Sample Demultiplexing (Parallelized):** The master variant list is then processed against each original sample to assign each mutation to its specific barcode of origin, generating the final output files for the visualization app.

This design, coupled with robust error handling and parallel execution, ensures reproducibility, scalability, and ease of use for large-scale genomic data analysis.

---

## 3. Detailed Script Analysis

### 3.1. Initial Setup and Configuration

The script begins with robust `set` commands for safe execution and an `error_handler` trap for detailed debugging. The most significant architectural change is the move to an external configuration model.

#### 3.1.1. Configuration Loading

The script no longer contains hard-coded paths or parameters. Instead, it loads them dynamically at runtime.

* **Argument Check:** The script first verifies that it was launched with exactly two arguments: the main config file and the barcode config file. This prevents execution with incomplete information.
* **Main Configuration (`pipeline.conf`):** It uses the `source` command to execute `pipeline.conf` in the current shell. This loads all system-specific variables, such as resource allocations (`jobs`, `threads_per_job`) and the path to the reference genome (`ref_fa`).
* **Barcode Configuration (`my_barcodes.cfg`):** The script reads the barcode file line by line, populating a Bash associative array named `barcodes`. This structure efficiently stores each barcode's name as a key and its sequence as the value. The loop includes logic to ignore empty lines or comments (lines starting with `#`), making the file flexible.

### 3.2. Workflow Stages (Functions)

The core logic is encapsulated within several Bash functions.

#### 3.2.1. Align & Filter (`align_and_filter`)

This function performs the initial alignment and filtering. The process involves converting the input BAM back to FASTQ, re-aligning with `bwa mem` using stringent mismatch and gap penalties, and using `tee` to split the output stream into two BAM files: one heavily filtered (`-F 4079 -q 60`) for consensus calling, and one unfiltered for mismatch extraction.

#### 3.2.2. BED-making & Base Counting (`make_bed_and_count`)

This function uses `bedtools bamtobed` to convert the filtered BAM into BED format. A crucial `awk` script then parses this BED file to find overlapping forward and reverse read pairs, defining the high-confidence consensus regions. The total length of these regions is summed and saved to a `.bases` file for the sample.

#### 3.2.3. Split & Extract Mismatches (`split_and_extract`)

This function extracts mismatch information from the unfiltered BAM. To enhance speed, it splits the BAM into several chunks. For each chunk, it uses `sam2tsv` to compare reads to the reference. An `awk` script then identifies and formats mismatches, creating a unique key for each event (`readname:chr:pos:ref:alt`).

#### 3.2.4. Join Strands & Finalize (`join_strands`)

This function is the core of the error-correction strategy. It aggregates all mismatch data for a single sample, sorts the forward and reverse strand mismatches separately, and uses the `join` command to find identical mismatch keys present in both files. This ensures that only mutations observed on both strands of the original DNA molecule are retained, effectively filtering out strand-specific sequencing errors.

#### 3.2.5. Demultiplex Sample (`demultiplex_sample`)

This is a new, critical function that performs the final demultiplexing stage for a single sample.

* **FASTA Conversion:** To improve efficiency, it converts the sample's BAM file to a FASTA format once at the beginning.
* **Barcode Searching:** It then loops through the `barcodes` array loaded from the config. For each barcode, it uses `seqkit grep` to search the FASTA file and generate a list of read names containing that barcode sequence.
* **Precise Filtering:** This list of read names is then used as a filter with `awk`. This is a significant improvement over simple `grep`, as `awk` is used to check for an exact match in the specific column containing the read name. This prevents false positives from partial matches in other columns.
    * It filters the master mutation list (`all.joined.2`) to create a barcode-specific `.muts` file.
    * It filters the sample's consensus BED file to create a barcode-specific `.bed` file.
* **Final Base Counting:** Finally, it calculates the total bases for the barcode-specific BED file, creating the final `.bases` output.

### 3.3. Main Execution Flow

The main execution block at the end of the script clearly orchestrates the three distinct parts of the pipeline.

* **Part 1: Core Mutation Pipeline:** Uses `GNU Parallel` to run the `process_sample` wrapper function for all input BAM files. This function executes steps 3.2.1 through 3.2.4 for each sample.
* **Part 2: Global Merging and Filtering:** After all parallel jobs are complete, this stage runs. It concatenates all per-sample results and applies a final, stringent set of `sort`, `uniq`, and `awk` filters based on strand orientation and mapping quality (MAPQ 60) for both reads supporting a variant. This produces the master mutation list, `v1-joined/all.joined.2`.
* **Part 3: Demultiplexing:** The script again uses `GNU Parallel` to run the `demultiplex_sample` function for all samples, efficiently distributing the master mutation list into the final barcode-specific output files located in the `indexing/` directory.
