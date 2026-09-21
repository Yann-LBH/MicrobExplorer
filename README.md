# 🧬 [MicrobExplorer] — Comprehensive Metagenomic & Functional Profiling Pipeline

![R](https://img.shields.io/badge/R->=%204.2.0-blue.svg)
![Phyloseq](https://img.shields.io/badge/Phyloseq-Integrated-green.svg)
![Shiny](https://img.shields.io/badge/Shiny-Interactive-orange.svg)
![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)

**MicrobExplorer** is end-to-end R-based pipeline designed to integrate, analyze, and visualize metagenomic output from **Kaiju**, **MEGAHIT**, and **Bakta**. It automates taxonomic and functional (KEGG) annotation matching, performs robust statistical workflows, and provides an interactive **Shiny** application for visual exploration and dynamic custom styling.

---

## 🔑 Key Features

- **Multi-Source Integration:** Seamlessly processes reads, contigs, and KEGG annotations from Kaiju, MEGAHIT, and Bakta outputs.
- **Automated Reference Matching:** Downloads and links up-to-date NCBI taxonomy and KEGG databases directly into `phyloseq` objects.
- **Statistical Framework:**
  - Multivariate analysis: PERMANOVA, Beta dispersion (`betadisper`), PCA, and CAP (Constrained Analysis of Principal Coordinates).
  - Differential abundance: Automated DESeq2 testing.
- **Rich Visualization Suite:** Generates ready-to-publish PCA plots, heatmaps, barplots, differential heatmaps, differential bar plots, CAP plots, and volcano plots.
- **Interactive Shiny App:**
  - Dynamic data filtering and real-time visualization.
  - Graphical customizer: Save custom plot aesthetics and export them directly back into the core pipeline.
- **Multilingual Support:** Dynamic titles and labels in English/French, easily extensible to other languages via a centralized JSON database.
- **Quality Control & Benchmarking:** Automated QC report generation and performance benchmarking.
- **Flexible Export:** Outputs publishable PDFs alongside lightweight `.parquet` and `.rds` files optimized for Shiny.

---

## 📁 Repository Architecture

```text
MicrobExplorer/
├── config/
│   ├── config.yaml          # Main pipeline configuration file
│   ├── samples.tsv          # Sample tracking mapping
│   └── i18n.json            # Dynamic titles & translation database
├── data/                    # Managed input directory
│   ├── raw/                 # Raw reads (kaiju), contigs, and kegg inputs
│   ├── metadata/            # Sample metadata (metadata.xlsx)
│   ├── taxonomy_ncbi/       # NCBI reference database (auto-downloaded)
│   ├── taxonomy_megahit/    # MEGAHIT taxNames conversions
│   ├── pathway_bakta/       # KEGG pathway reference (auto-downloaded)
│   └── physico_parameters/  # Optional physicochemical parameters
├── results/                 # Automated output directory
│   ├── treatment/           # Cleaned and normalized tables
│   ├── plots/               # Publication-ready PDF figures
│   ├── permanova/           # Statistical tables & outputs
│   ├── qc/                  # Quality control summaries
│   ├── rds/ & parquet/      # Binary exports for Shiny
│   └── audit/               # HTML audit report
├── run_pipeline.sh          # Shell wrapper to execute the core pipeline
├── run_shiny.sh             # Shell wrapper to launch the interactive GUI
└── README.md
```

---

## ⚙️ Configuration (`config/config.yaml`)

All pipeline behaviors, thresholds, and module activations are controlled via the `config/config.yaml` file. 

### Key Configuration Sections:

- **Language Support (`language`):** Set to `"en"` (English) or `"fr"` (French) for dynamic title generation across plots and summaries.
- **Module Toggles (`run_*`):** Seamlessly activate or deactivate specific pipeline components:
  - **Data Inputs:** `run_reads`, `run_contigs`, `run_kegg`, `run_physico`
  - **Analyses:** `run_permanova`, `run_deseq2`, `run_phyloseq`
  - **Visualizations:** `run_stackedbarplot_abundance`, `run_stackedbarplot_deseq2`, `run_heatmap_abundance`, `run_heatmap_deseq2`, `run_pca`, `run_volcano`
  - **Utilities:** `run_taxonomy`, `run_qc`, `run_shiny`, `run_benchmarks`
- **Filtering & Quality Thresholds:**
  - **Reads:** Set minimum raw count limits (`count_threshold`) and exclude uninformative taxa or contaminants via a custom `noise` list.
  - **Contigs:** Filter assembly metrics using `length_threshold`, `abundance_threshold`, and `rpkm_threshold`.
  - **KEGG Functional Analysis:** Adjust statistical cutoffs using `pvalue_threshold` and `lfc_threshold`.
- **Statistical Parameters:**
  - **PERMANOVA:** Define experimental factors (`effect`), distance metrics (`distance_method`), CLR transformation (`use_clr`), and permutations (`permutation`).
  - **DESeq2:** Configure experimental design contrasts (`contrast`, `ref`), test type (`Wald` or `LRT`), and fit type.
- **Plotting Aesthetics & Export (`plots`):**
  - Custom global dimensions for print-ready PDFs (`pdf_size`).
  - Fine-tune color palettes (`palette: "turbo"`), themes (`theme_minimal`), font sizes, and top-$N$ feature display limits.

---

## 📋 Input Metadata & Sample Files

To map raw outputs to experimental variables, MicrobExplorer requires two configuration files in the `data/` and `config/` directories:

### 1. Sample Mapping (`config/samples.tsv`)
A file naming each sample ID you want to use:

```text
sample
T1_ctrl_261501
T1_261501
```

### 2. Experimental Metadata (`data/metadata/metadata.xlsx`)
An Excel sheet containing experimental variables used for PERMANOVA, DESeq2 contrasts, and plot grouping (`sample_id`, `group`, `name`, `date`):

```text
| sample_id | group | name | date |
| :--- | :--- | :--- | :--- | :--- |
| **T1_ctrl_261501** | Ctrl | T1_Ctrl | 2026/01/15 |
| **T1_261501** | T1 | T1 | 2026/01/15 |
```

### 3. Physicochemical Metadata (`data/physico_parameters/parameters.xlsx`)
An Excel sheet containing physicochemical variables used for ACP, and plot grouping (`date`, `name`, `ph`, etc.):

```text
| date | name | pH | AGV mg/l |
| :--- | :--- | :--- | :--- | :--- |
| **2026/01/15** | T1_Ctrl | 6.8 | 42 |
| **2026/01/15** | T1 | TD1_T1 | 7.2 | 34 |
```

---

## 🛠️ Environments & Dependencies

MicrobExplorer relies on two dedicated Conda environments (`r-env` and `py-env`) to ensure full reproducibility and isolate dependencies.

### Prerequisites
- Unix-based OS (Linux or macOS)
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- `git`

### 1. Setup Environments

Clone the repository and build the Conda environments:

```bash
# Clone the repository
git clone [https://github.com/Yann-LBH/MicrobExplorer.git](https://github.com/Yann-LBH/MicrobExplorer.git)
cd MicrobExplorer

# Create the Conda environments from the environment files
conda env create -f config/r-env.yml
conda env create -f config/py-env.yml

# Make execution scripts executable
chmod +x run_pipeline.sh run_shiny.sh

# Rendre les scripts shell exécutables
./run_pipeline.sh --config config/config.json

# Shiny app
./run_shiny.sh