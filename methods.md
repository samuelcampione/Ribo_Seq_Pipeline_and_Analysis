# __Method__
## Processing Raw Data ([Shell script pipeline](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/ribo_profiling_pipeline_stdin.sh))
- Acquire raw reads data from GEO
- Trim adapters and filter low quality reads (cutadapt)
- Deduplicate (seqkit rmdup)
- Remove unique barcodes (cutadapt)
- Filter rRNA contamination (bowtie2)
- Align and quantify reads (STAR)
## Statistical Analysis ([R script for DE analysis](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/GIR2%20vs%20WT%20differential%20analysis.R))
- Filter Out Lowly Expressed Genes
- Normalize counts and correct for batch effects (edgeR and ComBat-seq)
- Differential Expression Analysis (edgeR)
- Identify differentially expressed genes and extract log2 fold changes of specific genes

## Dependencies
See [environment.yml](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/environment.yml) for package versions.
