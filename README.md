<p align="center">
  <img src="assets/PHASE_Suite_logo.png" alt="PHASE Suite Logo" width="200"/>
</p>

<h1 align="center">PHASE Suite</h1>

<p align="center">
  A command‑line toolkit for somatic mutation detection from error-corrected next-generation sequencing PacBio HiFi data in genetic toxicology.
</p>

<p align="center">
  <a href="https://github.com/PHASE-Suite/PHASE-Suite/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License"></a>
  <a href="https://github.com/PHASE-Suite/PHASE-Suite/releases"><img src="https://img.shields.io/github/v/release/PHASE-Suite/PHASE-Suite?label=Latest%20Release" alt="Latest Release"></a>
  <a href="https://doi.org/10.5281/zenodo.16624362"><img src="https://zenodo.org/badge/1028706985.svg" alt="DOI"></a>
  <a href="https://phase-suite.readthedocs.io/"><img src="https://img.shields.io/readthedocs/phase-suite" alt="Documentation Status"></a>
</p>

---

## 📖 Introduction

Welcome to the **PHASE Suite**! **P**acBio **H**iFi **A**nalysis of **S**omatic **E**vents.

In disciplines such as genetic toxicology, accurately identifying DNA mutations is crucial for determining the genotoxic potential of chemicals and evaluating cancer risk. Although next-generation sequencing (NGS) has revolutionized genetic research, its effectiveness is limited by a major flaw: a built-in error rate ranging from 0.1% to 1% per base. This level of technical noise often overwhelms the actual mutation frequency, which may be as low as one in 100 million (1×10⁻⁸). Consequently, rare but important mutations relevant to cancer studies and safety testing are often obscured by sequencing errors and go undetected.

To address this issue, error-corrected NGS (ecNGS) technologies have become an essential advancement. These methods—including PacBio HiFi sequencing—enhance accuracy by sequencing both strands of DNA redundantly, significantly reducing background errors and allowing for the reliable detection of extremely rare mutations.

This project provides a powerful and accessible set of tools that implements this ecNGS strategy for PacBio HiFi bam files. By computationally pairing the forward and reverse strands of a single DNA molecule, PHASE filters out random errors to confidently call true mutations.

The suite currently includes:

**PHASE Pipeline** ([`phase-script.sh`](phase-script.sh)): A robust bash script that automates the entire workflow, from raw PacBio BAM files to a final, high-confidence list of somatic mutations.

**Visualization Script** (Coming Soon): A Python tool for generating publication-quality plots of mutational signatures from the pipeline's output.

---

## ✨ Features

* **High-Fidelity Mutation Calling**: Achieves high sensitivity and precision by requiring that a mutation be present on both the forward and reverse strands of a single DNA molecule.
* **Automated Workflow**: A single script orchestrates a multi-step process involving alignment, filtering, strand-pairing, and variant calling (**mismatches**).
* **Parallel Processing**: Built with GNU Parallel to efficiently process multiple samples across many CPU cores, significantly reducing runtime.
* **Comprehensive Filtering**: Implements several layers of filtering based on **mapping quality**, **strand agreement**, and the removal of **non-unique mutation events** to ensure high-confidence results.
* **Clear Outputs**: Generates well-structured output files, including a final tab-separated file of mutations ready for downstream analysis.
* **(Coming Soon) Mutational Signature Analysis**: The accompanying Python script will allow for easy visualization of the mutational landscape (e.g. SBS signatures).

---

## 🔬 In-Depth Technical Analysis

For a comprehensive, command-by-command breakdown of the pipeline's methodology and scientific rationale, please see our detailed technical analysis.

[**Read the Full Technical Analysis Here**](TECHNICAL_ANALYSIS.md)

---

## 🚀 Installation & Quick Start

The recommended way to install PHASE Suite and its dependencies is by using the provided Conda environment file.

**1. Clone the Repository**
```bash
git clone [https://github.com/PHASE-Suite/PHASE-Suite.git](https://github.com/PHASE-Suite/PHASE-Suite.git)
cd PHASE-Suite
```

**2. Create and Activate the Conda Environment**

This command uses the `environment.yml` file to automatically install all required software.
```bash
conda env create -f environment.yml
conda activate phase_suite_env
```

**3. Configure Pipeline Paths**

Before running, you must open `phase-script.sh` in a text editor and provide the correct paths for your reference genome and the `sam2tsv.jar` file.
```bash
# Open phase-script.sh and edit these lines near the top:
ref_fa="/path/to/your/reference_genome.fa"
sam2tsv_jar="/path/to/your/jvarkit/sam2tsv.jar"
```

**4. Run the Pipeline**

Place your input BAM files in the `original-bams/` directory and execute the script:
```bash
bash phase-script.sh
```

---

## 📂 Repository Structure

An overview of the key files and directories in this project:

```
PHASE-Suite/
├── phase-script.sh       # <== The main executable bash pipeline
├── TECHNICAL_ANALYSIS.md # <== In-depth explanation of the script's logic
├── environment.yml       # <== Conda environment file for dependencies
├── original-bams/        # <== Directory for your input BAM files
└── README.md             # <== You are here!
```
