# mutSOMA: Somatic Mutation Rate Estimation in Trees

## Project Structure

- `R_Wrapper.R` – Main function to execute the pipeline
- `makeVCFpedigreeTEMPv2.R` – Prepares pedigree data from VCF files
- `makePHYLO.R` - Calculates sample divergence times
- `bootSOMA.R` – Bootstrapping for confidence intervals
- `additional_files/` – Sample metadata and reference genome FASTA file
- `pedigree_files/` - Saves pedigree data
- `out/`- Saves estimation results

---

## Requirements

Install the following R packages:

```r
install.packages("vcfR")
install.packages("expm")
install.packages("reticulate")

# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

The reference FASTA file is compressed as a `.7z` archive. Please extract it before using it. 

#### On Linux/macOS:
```bash
7z x PtrichocarpaStettler14_532_v1.0.fa.7z
