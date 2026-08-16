# Equine Uterine Microbiome — QIIME2 Workflow

### DNA- and RNA-based 16S rRNA analysis of the low-biomass equine uterine microbiome

This repository contains the **QIIME2 workflow** developed to compare DNA- and RNA-based 16S rRNA V3–V4 amplicon sequencing approaches for the characterization of the low-biomass bacterial microbiome of the equine uterus.

The workflow accompanies the study:

**Dyroff, A.I., López-Valiñas, Á., Magalhaes, H.B. et al. (2025).**  
*Comparison of RNA- and DNA-based 16S amplicon sequencing to find the optimal approach for the analysis of the uterine microbiome.*  
Scientific Reports 15, 17037.  
https://doi.org/10.1038/s41598-025-00969-5

> **Important:** The additional files required to reproduce the analysis — including the trained taxonomic classifier, SILVA reference database, manifest, metadata and archived workflow resources — are available in the associated **Zenodo repository:**  
> https://doi.org/10.5281/zenodo.13847011

---

## Workflow overview

```text
Illumina paired-end reads
          │
          ▼
      QIIME2 import
          │
          ▼
 Cutadapt primer removal
          │
          ▼
     DADA2 denoising
          │
          ▼
 Feature table + ASVs
          │
          ▼
 Non-prokaryotic sequence removal
          │
          ▼
 Taxonomic classification
    SILVA / scikit-learn
          │
          ▼
 Chloroplast + mitochondrial removal
          │
          ▼
     ASV filtering
  abundance + prevalence
          │
          ├────────────────────┐
          ▼                    ▼
   Taxonomic tables      MAFFT alignment
   relative/absolute            │
                               ▼
                         Sequence masking
                               │
                               ▼
                         FastTree phylogeny
                               │
                               ▼
                       Midpoint-rooted tree
                               │
                               ▼
                       Diversity analyses
                     ┌─────────┴─────────┐
                     ▼                   ▼
              Alpha diversity       Beta diversity
          Shannon / Simpson / Chao1     Jaccard
                                         │
                                         ▼
                                 PERMANOVA / ADONIS
```

The workflow additionally generates taxonomic abundance tables across multiple taxonomic levels for downstream analyses.

---

## Analysis strategy

Raw paired-end Illumina reads are imported into QIIME2 and primers are removed with **Cutadapt** before denoising and ASV inference with **DADA2**.

Representative sequences are filtered against the SILVA reference to remove non-prokaryotic sequences and subsequently classified using a trained **SILVA / scikit-learn classifier**. Mitochondrial and chloroplast sequences are excluded before downstream analysis.

Low-abundance ASVs are filtered using the strategy implemented in the original study: ASVs supported by at least **5 reads in at least 5 samples** are retained, while otherwise excluded ASVs can be recovered when supported by at least **20 reads in one sample**.

The filtered ASV dataset is subsequently used for phylogenetic reconstruction and diversity analyses. Representative sequences are aligned with **MAFFT**, hypervariable alignment positions are masked, and a phylogenetic tree is reconstructed using **FastTree** and midpoint rooting.

Alpha diversity analyses include **Shannon, Simpson and Chao1** metrics. Beta diversity analysis includes **Jaccard distances**, group comparisons and **PERMANOVA / ADONIS**.

---

## Software

The original analysis was performed using:

- **QIIME2 2024.5**
- macOS Ventura 13.4
- Apple M1 Pro
- `osx-arm64`

Additional tools used through QIIME2 or the command line include **DADA2, Cutadapt, VSEARCH, MAFFT, FastTree and BIOM**.

---

## Data and resources

Raw sequencing data are publicly available through NCBI:

**BioProject:** [PRJNA1165023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1165023)

Additional resources required to reproduce the workflow are available on Zenodo:

**Zenodo:** [10.5281/zenodo.13847011](https://doi.org/10.5281/zenodo.13847011)

These resources include the trained classifier, SILVA reference sequences, manifest, metadata and archived analysis files and should be downloaded before running the workflow.

---

## Running the workflow

1. Install **QIIME2 2024.5**.
2. Download the FASTQ files from **[BioProject PRJNA1165023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1165023)**.
3. Download the required reference and metadata files from the associated **[Zenodo repository](https://doi.org/10.5281/zenodo.13847011)**.
4. Update the FASTQ paths in `manifest.csv` for your local system.
5. Run:

```bash
bash pipeline.sh
```

The complete analysis is provided in `pipeline.sh`, following the same sequence of processing, filtering, phylogenetic and diversity-analysis steps used in the associated study.

---

## Citation

If you use this workflow, please cite the associated publication:

> Dyroff, A.I., López-Valiñas, Á., Magalhaes, H.B. et al. (2025). Comparison of RNA- and DNA-based 16S amplicon sequencing to find the optimal approach for the analysis of the uterine microbiome. *Scientific Reports* 15, 17037. https://doi.org/10.1038/s41598-025-00969-5

The archived workflow and associated analysis resources are available from Zenodo:

> López-Valiñas, Á. (2024). *Unveiling the optimal approach for the analysis of the low biomass bacterial microbiome in the equine uterus: Comparison of RNA- versus DNA-based 16S rRNA V3-4 amplicon*. Zenodo. https://doi.org/10.5281/zenodo.13847011

---

## Author

**Álvaro López-Valiñas, PhD**  
Bioinformatics · Microbiome · Genomics  
ORCID: 0000-0002-7492-5108
