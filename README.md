# hnRNP R knockout RNA-seq analysis

Analysis code for:

> **Loss of hnRNP R Disrupts SMN-Associated Ribonucleoprotein Homeostasis**
> Changhe Ji. *Biology* (MDPI), manuscript biology-4683043.

This repository holds every script used to analyse the RNA-sequencing data behind
Figures 5, 6 and 7 and the corresponding supplementary figures. Nothing here is
re-implemented or tidied for presentation: these are the scripts as they were run,
so that the published figures can be traced to the code that produced them.

## Data

| What | Where |
|---|---|
| Raw and processed sequencing data | GEO accession **GSE202720** (Ji *et al.*, *EMBO Rep* 2022) |
| Analysis output behind every plotted panel | `Source_Data_biology-4683043.xlsx`, published with the paper |
| Cell lines | Prime-edited HNRNPR<sup>+/+</sup> and HNRNPR<sup>−/−</sup> HeLa clones from Ji *et al.*, *EMBO Rep* 2022 |

The sequencing data are **not** redistributed here. Download them from GEO; `data/`
is where the scripts expect them.

## Pipeline

Upstream processing (not in this repository, run from the command line):

| Step | Tool | Version |
|---|---|---|
| Adapter and quality trimming | Cutadapt | 2.5 |
| Alignment | STAR | 2.7.2b |
| Gene-level counting | featureCounts (Subread) | 1.6.4 |
| Genome / annotation | GRCh38.p13, RefSeq 109.20210226 | |
| Alternative splicing | Shiba | |
| Alternative polyadenylation | QAPA (GENCODE v50) | |
| Transcription start site usage | DATTSS | see `shell/12_TSS_DATTSS_commands.sh` |

Downstream analysis in this repository:

| Script | What it does | Figures |
|---|---|---|
| `R/01_differential_expression_DESeq2.R` | DESeq2 differential expression, volcano plot | 5B |
| `R/02_CPM_TPM_tables.R` | CPM and TPM tables | — |
| `R/03_heatmap_VST.R` | VST-transformed expression heatmap | 5C |
| `R/04_heatmap_top30_up_down.R` | Top-30 up / down heatmap | 5C |
| `R/05_GO_ORA_GSEA_clusterProfiler.R` | GO over-representation **and** GSEA, clusterProfiler, expressed-gene background, Benjamini–Hochberg | **5D** |
| `R/06_GO_results_tables.R` | GO result tables | 5D |
| `R/07_GO_next_level_outputs.R` | Additional GO outputs | supplementary |
| `R/08_GO_next_level_revised.R` | Revised version of the above | supplementary |
| `R/09_enrichment_gprofiler2.R` | gprofiler2 `gost` enrichment, cross-check | supplementary |
| `R/10_splicing_PSI_PCA.R` | PCA of percent-spliced-in values | 6A |
| `R/11_APA_QAPA.R` | QAPA poly(A) site usage, 3′UTR length, volcano and heatmap | 7A–C |
| `shell/12_TSS_DATTSS_commands.sh` | DATTSS command lines for tandem TSS usage | 7D |

## Which script produced Figure 5D

Figure 5D is produced by **`R/05_GO_ORA_GSEA_clusterProfiler.R`**. This matters because an
earlier draft of the manuscript showed a GO analysis run without a proper background gene
set, which returned adhesion, MAPK and extracellular-matrix terms and disagreed with the
GSEA of the same data. Run with the expressed genes as background, over-representation
analysis and GSEA agree: ribonucleoprotein complex biogenesis, ribosome biogenesis, rRNA
and mRNA processing and RNA splicing, with `cytoplasmic stress granule` among the Cellular
Component terms. Figure 5D in the published paper is the latter analysis.

## Known issues in the code, left as they are

These are recorded rather than silently corrected, so that the archive matches what was run:

- **Absolute paths.** Every script carries the author's own Windows paths
  (`H:/JI3/HNRNPR/workdir_v8.1.0/...`) or WSL paths (`/mnt/h/JI3/HNRNPR/...`). To run the
  code elsewhere, point them at a local directory, for example:
  `sed -i 's|H:/JI3/HNRNPR/workdir_v8.1.0|/your/path|g' R/*.R`
- **`R/02_CPM_TPM_tables.R`** contains a GO block that reads `H:/7SK/7SK/workdir1/results/expression/counts.txt`
  — a counts file from a different project, left over from an earlier script this one was
  adapted from. That block selected terms by the top 30 genes ranked on mean difference,
  with no statistical test and no background set. **No figure in this manuscript uses its
  output.** All GO results in the paper come from `R/05_GO_ORA_GSEA_clusterProfiler.R`.
- **No `renv.lock`.** Package versions were not pinned at the time of the analysis. The
  packages used are listed below; the versions are those current in the author's R
  installation when the scripts were run.

## Requirements

R, with:

```
DESeq2  clusterProfiler  org.Hs.eg.db  enrichplot  ReactomePA  gprofiler2
ggplot2  ggrepel  patchwork  aplot  scales  pheatmap
dplyr  tidyr  tibble  readr  readxl  openxlsx  data.table  stringr  matrixStats
```

```r
install.packages(c("ggplot2","ggrepel","patchwork","aplot","scales","pheatmap",
                   "dplyr","tidyr","tibble","readr","readxl","openxlsx",
                   "data.table","stringr","matrixStats","gprofiler2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","clusterProfiler","org.Hs.eg.db","enrichplot","ReactomePA"))
```

## Citing

If you use this code, please cite the paper and this archive. The archived release has its
own DOI; see `CITATION.cff`.

## Licence

MIT — see `LICENSE`.
