RNAseq 6w vs 12w +/- MMTV-MYC mouse merged data June-September 2024

# First attempt
k693_rerun_analysis_GS.R --> analysis with interaction myc-status vs timepoint

  - coldata setup: 

from: data/FULL.DAT.COL.DATA.txt
columns:
- sample
- group
- myc_status
- timepoint
- prefix
- suffix
- numeric_part
- alphabetical_part

- counts: cts
- issue?: 12W has about half total counts comapred to 6W

- dds: design = ~ timepoint * myc_status) #this expands to ~ timepoint + myc_status + timepoint:myc_status
  
- coverage: 19K genes >10 counts
  
- issue?: data are noisy - none of the variance stabilising transformations, followed by distance measures or PCA ("log2(x + 1)", "vst", "rlog") reveal grouping by the actual experimental groups
PC1 and 2 top genes hint to immune cells...
  
- further exploration with variancePartition: only small part of the variation is explained by any of the parameters: see coldata: (1|myc_status) + (1|timepoint) + (1|numeric_part) + (1|alphabetical_part) + (1|myc_status) + (1|timepoint) (1|prefix) (1|suffix)

- moving on for DGE anyway: 

- altogether this analysis provides a list, using the interaction terms to create lists of genes whihc are significant at 6W and 12W. The genes are used in Gprofiler, but it is is difficult to interpret the changes.

- I have also tried a Cytoscpae analysis to visualise the gene lists, but it is not yet completed and/or informative.

- Conclusion: this file is used as the base of the second analysis where both the interactiion and group-based design is used to create gene lists, to understand the difference in the Myc effect between 6W and 12W timepoints.

# Second attempt
Myc_timecourse_analysis_GS.R


  

