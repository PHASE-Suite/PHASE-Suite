<h1 align="center">PHASE Suite (PacBio HiFi Analysis of Somatic Events)</h1>

<p align="center">
A complete toolkit for high-fidelity somatic mutation detection, analysis, and interactive visualization from PacBio HiFi data.
</p>

<p align="center">
<a href="https://github.com/PHASE-Suite/PHASE-Suite/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License"></a>
<a href="https://github.com/PHASE-Suite/PHASE-Suite/releases"><img src="https://img.shields.io/github/v/release/PHASE-Suite/PHASE-Suite?label=Latest%20Release" alt="Latest Release"></a>
<a href="https://doi.org/10.5281/zenodo.16624362"><img src="https://zenodo.org/badge/1028706985.svg" alt="DOI"></a>
</p>

---

## 📖 Introduction

Welcome to the **PHASE Suite**!

Accurate detection of DNA mutations is essential in fields like genetic toxicology, where it's used to assess the genotoxicity of chemicals and evaluate cancer risk. While next-generation sequencing (NGS) has transformed genetic research, its utility is limited by an inherent error rate of 0.1% to 1% per base. This background noise often exceeds the true mutation rate, which can be as low as 1 in 100 million (1×10⁻⁸), making it difficult to detect rare but biologically important mutations.

Error-corrected NGS (ecNGS) technologies, such as **PacBio HiFi sequencing**, overcome this limitation by redundantly sequencing both DNA strands, dramatically reducing errors and enabling the reliable detection of ultra-rare mutations.

**PHASE Suite** is a complete, user-friendly bioinformatics toolkit for analyzing PacBio HiFi data. By computationally pairing the forward and reverse strands of individual DNA molecules, it filters out random sequencing errors, enabling accurate identification and visualization of true low-frequency mutations.

🏛️ **Suite Components**
The suite includes two main components:

**PHASE Pipeline** ([`phase_pipeline.sh`](phase_pipeline.sh)): A robust Bash script that automates the entire workflow, from raw PacBio BAM files to a final, demultiplexed list of high-confidence somatic mutations.

**PHASE Explorer** ([`phase_explorer.py`](phase_explorer.py)): An interactive Streamlit application for generating publication-quality plots, performing statistical analysis, and exploring mutational signatures from the pipeline's output.

---

## ✨ Features

* **High-Fidelity Mutation Calling**: Achieves high sensitivity and precision by requiring that a mutation be present on both the forward and reverse strands of a single DNA molecule.
* **Complete End-to-End Workflow**: A unified suite that takes you from raw BAM files to interactive, publication-ready figures.
* **Modular Configuration**: Pipeline settings (paths, threads) and experimental parameters (barcodes) are managed in simple text files, separate from the core code.
* **Parallel Processing**: Built with GNU Parallel to efficiently process multiple samples across many CPU cores, significantly reducing runtime.
* **Interactive Visualization**: The PHASE Explorer app provides dynamic plots with hover-to-view details, zooming, and advanced statistical analysis (ANOVA, t-tests).
* **Comprehensive Filtering**: Implements several layers of filtering based on mapping quality, strand agreement, and uniqueness to ensure high-confidence results.

---

## 🎓 Publications & Citations

The PHASE Suite pipeline and its underlying methodology have been successfully used in the following peer-reviewed studies:

* **Mutation accumulation following extended exposure of human HepaRG cells to a genotoxic carcinogen**. Journal of Environmental Science and Health, Part C. Published: 2025-08-22. [DOI: 10.1080/26896583.2025.2548146](https://doi.org/10.1080/26896583.2025.2548146)

* **Mutation Frequencies in TK6 and L5178Y Cells: Implications for Error-Corrected Sequencing**. Environmental and molecular mutagenesis. Published: 2025-07-15. [DOI: 10.1002/em.70024](https://doi.org/10.1002/em.70024)

* **Assessment of in vivo chemical mutagenesis by long-read sequencing**. Toxicological Sciences. Published: 2024-11-01. [DOI: 10.1093/toxsci/kfae104](https://doi.org/10.1093/toxsci/kfae104)

* **Evaluating the mutagenicity of N-nitrosodimethylamine in 2D and 3D HepaRG cell cultures using error-corrected next generation sequencing**. Archives of Toxicology. Published: 2024-04-07. [DOI: 10.1007/s00204-024-03731-4](https://doi.org/10.1007/s00204-024-03731-4)

* **Unbiased whole genome detection of ultrarare off-target mutations in genome-edited cell populations by HiFi sequencing**. Environmental and molecular mutagenesis. Published: 2023-07-24. [DOI: 10.1002/em.22566](https://doi.org/10.1002/em.22566)

* **Evaluation of the mutagenic effects of Molnupiravir and N4-hydroxycytidine in bacterial and mammalian cells by HiFi sequencing**. Environmental and molecular mutagenesis. Published: 2022-08-01. [DOI: 10.1002/em.22510](https://doi.org/10.1002/em.22510)

* **Genome-wide detection of ultralow-frequency substitution mutations in cultures of mouse lymphoma L5178Y cells and Caenorhabditis elegans worms by PacBio sequencing**. Environmental and molecular mutagenesis. Published: 2022-02-01. [DOI: 10.1002/em.22473](https://doi.org/10.1002/em.22473)

* **PacBio sequencing detects genome-wide ultra-low-frequency substitution mutations resulting from exposure to chemical mutagens**. Environmental and molecular mutagenesis. Published: 2021-09-02. [DOI: 10.1002/em.22462](https://doi.org/10.1002/em.22462)

    
---

## 🔬 In-Depth Technical Analysis

For a comprehensive, command-by-command breakdown of the methodology and scientific rationale, please see our detailed technical analysis documents:

[**Analysis of the Data Processing Pipeline**](TECHNICAL_ANALYSIS.md)

[**Analysis of the Visualization App**](TECHNICAL_ANALYSIS_APP.md)

---

## 🚀 Installation & Quick Start

The recommended way to install the PHASE Suite and all its dependencies is by using the provided Conda environment file.

**1. Clone the Repository**
```bash
git clone https://github.com/PHASE-Suite/PHASE-Suite.git
cd PHASE-Suite
```

**2. Create and Activate the Conda Environment**
This single command uses the environment.yml file to install all required software for both the pipeline and the visualization app.
```bash
conda env create -f environment.yml
conda activate phase_suite_env
```

**3. Prepare Your Data and Configurations**

* **Place** your input `.bam` files into the `original-bams/` directory.

   * **Important**: Name your files using the format `Identifier_Condition.bam` (e.g., `TK6_10uM-N2.bam`) for the visualization app to correctly parse them.

* **Edit** `pipeline.conf` and set the `ref_fa` variable to the absolute path of your BWA-indexed reference genome.

* **Edit or create** `my_barcodes.cfg` to list the specific barcodes and sequences for your experiment.

**4. Run the Pipeline**
Execute the pipeline script, providing your two configuration files as arguments.
```bash
bash phase_pipeline.sh pipeline.conf my_barcodes.cfg
```

**5. Launch the Visualization App**
Once the pipeline is finished, launch the PHASE Explorer app with Streamlit.

```
streamlit run phase_explorer.py
```

The app will open in your web browser. Follow the on-screen instructions to point it to your indexing/ folder and reference genome to start visualizing your results.

---

## 📂 Repository Structure

An overview of the key files and directories in this project:

```
PHASE-Suite/
├── phase_pipeline.sh           # <== The main executable bash pipeline
├── phase_explorer.py           # <== The interactive Streamlit visualization app
├── pipeline.conf               # <== Main configuration for the pipeline
├── my_barcodes.cfg             # <== Example barcode file for demultiplexing
├── environment.yml             # <== Conda file to install all dependencies
├── original-bams/              # <== Directory for your input BAM files
├── TECHNICAL_ANALYSIS.md       # <== In-depth explanation of the pipeline's logic
├── TECHNICAL_ANALYSIS_APP.md   # <== In-depth explanation of the app's logic
└── README.md                   # <== You are here!
```
