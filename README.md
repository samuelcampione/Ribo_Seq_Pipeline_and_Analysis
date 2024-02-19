# Bioinformatics Analysis of Mined Ribo-seq Data

This is a bioinformatics workflow for analysis of raw ribosome profiling data. 
To ensure my analysis is sound, I reproduced key results of Egorov et al. (2021)'s publication using their available Ribo-seq data.





__Samples Used__
- GIR2Δ Strain
- WT Strain
- PUB1Δ Strain (for batch correction)



__Method__
- Processing Raw Data ([Shell script pipeline](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/ribo_profiling_pipeline_stdin.sh))
  - Acquire raw reads data from GEO
  - Trim adapters and filter low quality reads (cutadapt)
  - Deduplicate (seqkit rmdup)
  - Remove unique barcodes (cutadapt)
  - Filter rRNA contamination (bowtie2)
  - Align and quantify reads (STAR)
- Statistical Analysis ([R script for DE analysis](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/GIR2%20vs%20WT%20differential%20analysis.R))
  - Filter Out Lowly Expressed Genes
  - Normalize counts and correct for batch effects (edgeR and ComBat-seq)
  - Differential Expression Analysis (edgeR)
  - Identify differentially expressed genes and extract log2 fold changes of specific genes




Data is available at [GSE185458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185458) and [GSE185286](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185286).


Egorov, A. A., Makeeva, D. S., Makarova, N. E., Bykov, D. A., Hrytseniuk, Y. S., Mitkevich, O. V., Urakov, V. N., Alexandrov, A. I., Kulakovskiy, I. V., & Dmitriev, S. E. (2021). Ribo-Seq and RNA-Seq of TMA46 (DFRP1) and GIR2 (DFRP2) knockout yeast strains. F1000Research, 10, 1162. https://doi.org/10.12688/f1000research.74727.1
