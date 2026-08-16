# Equine Uterine Microbiome — QIIME2 Workflow

Bioinformatics workflow for the analysis of **DNA- and RNA-based 16S rRNA V3–V4 amplicon sequencing data** from the low-biomass bacterial microbiome of the equine uterus.

This repository contains the QIIME2 workflow developed for the study:

**Unveiling the optimal approach for the analysis of the low biomass bacterial microbiome in the equine uterus: Comparison of RNA- versus DNA-based 16S rRNA V3-4 amplicon**

## Workflow overview

The pipeline performs:

* Import of paired-end Illumina sequencing reads
* Primer and adapter removal using **Cutadapt**
* ASV inference and denoising using **DADA2**
* Removal of non-prokaryotic sequences
* Taxonomic classification using a trained **SILVA / scikit-learn classifier**
* Removal of mitochondrial and chloroplast sequences
* Abundance-based ASV filtering
* Multiple sequence alignment using **MAFFT**
* Phylogenetic reconstruction using **FastTree**
* Alpha diversity analysis
* Beta diversity analysis
* PERMANOVA / ADONIS analysis
* Generation of taxonomic abundance tables at multiple taxonomic levels

## Software

The original analysis was performed using:

* **QIIME2 2024.5**
* macOS Ventura 13.4
* Apple M1 Pro
* `osx-arm64`

Additional tools are called through QIIME2 and the command line, including DADA2, Cutadapt, VSEARCH, MAFFT, FastTree and BIOM.

## Data availability

Raw sequencing data are publicly available through NCBI:

**BioProject:** [PRJNA1165023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1165023)

The trained classifier, SILVA reference database, manifest, metadata and archived version of this workflow are available on Zenodo:

**Zenodo DOI:** [10.5281/zenodo.13847011](https://doi.org/10.5281/zenodo.13847011)

## Running the workflow

1. Install **QIIME2 version 2024.5**.

2. Download the FASTQ files from BioProject PRJNA1165023.

3. Download the required reference files from the associated Zenodo repository.

4. Update the paths in `manifest.csv` to point to the FASTQ files on your system.

5. Run:

```bash
bash pipeline.sh
```

## Repository contents

```text
equine-uterine-microbiome/
├── README.md
└── pipeline.sh
```

Reference databases, classifiers and additional project files are hosted on Zenodo rather than duplicated in this repository.

## Citation

If you use this workflow, please cite the associated Zenodo record:

> López-Valiñas, Á. (2024). *Unveiling the optimal approach for the analysis of the low biomass bacterial microbiome in the equine uterus: Comparison of RNA- versus DNA-based 16S rRNA V3-4 amplicon*. Zenodo. https://zenodo.org/records/13847011

## Author

**Álvaro López-Valiñas**
ORCID: [0000-0002-7492-5108](https://orcid.org/0000-0002-7492-5108)

