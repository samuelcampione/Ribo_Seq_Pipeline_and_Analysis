# Bioinformatics Analysis of Mined Ribo-seq Data

## Overview
This is a bioinformatics workflow for analysis of raw ribosome profiling data. 
To ensure my analysis is sound, I reproduced key results of Egorov et al. (2021)'s publication using their available Ribo-seq data.

__Samples Used:__
- GIR2Δ Strain
- WT Strain
- PUB1Δ Strain (for batch correction)

__Genes checked:__ 
- STR3, SNF3 (Upregulated in Egorov et al. (2021))
- YNL011C, snR86 (Downregulated in Egorov et al. (2021))
- SPO24 (Unchanged in Egorov et al. (2021))




## __Method__ 
-  Processing Raw Data ([pipeline.sh](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/pipeline.sh))
-  Statistical Analysis ([analysis.R](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/analysis.R))
- See [methods.md](https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/methods.md) for additional detailed steps



## Results

<table>
  <tr>
    <td>
      <img src="https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/visualizations/figure1B.png" width="400"/>
      <br>
      Results from DE analysis.
    </td>
    <td>
      <img src="https://github.com/samuelcampione/Ribo_Seq_Pipeline_and_Analysis/blob/main/visualizations/gir2_vs_wt_logFC.png" width="400"/>
      <br>
      Figure 1B from Egorov et al. (2021).
    </td>
  </tr>
</table>

- Upregulated similar to Egorov et al.
  - STR3: 	   1.6 logFC
  - SNF3: 	   1.4 logFC
- Downregulated similar to Egorov et al.
  - snR86:	  -2.1 logFC
  - YNL011C:  -1.8 logFC
- SPO24 is not significantly changed.

My results concur with the findings in the publication (Egorov et al. (2021)).



## Data Availability
Data is available at [GSE185458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185458) and [GSE185286](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185286).


_Egorov, A. A., Makeeva, D. S., Makarova, N. E., Bykov, D. A., Hrytseniuk, Y. S., Mitkevich, O. V., Urakov, V. N., Alexandrov, A. I., Kulakovskiy, I. V., & Dmitriev, S. E. (2021). Ribo-Seq and RNA-Seq of TMA46 (DFRP1) and GIR2 (DFRP2) knockout yeast strains. F1000Research, 10, 1162. https://doi.org/10.12688/f1000research.74727.1_
