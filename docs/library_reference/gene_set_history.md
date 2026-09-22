**Myc-mouse/main**
1. in scripts/archive_main_pipeline/Myc_timecourse_analysis_GS_sandbox.R: 
- gene_sets_list object created (l612), this has:
    - custom selected mitocarta pathways (/Users/gs/G/data/MK_myc_2022/myc_mouse/data/mitocarta_pathways.csv)
    - myc signature compendium downloaded from Felsher paper: /Users/gs/G/data/MK_myc_2022/myc_mouse/data/myc_signature_genesets.gmx
    - Felsher integrative gene set added (/Users/gs/G/data/MK_myc_2022/myc_mouse/data/felsher_integrative_signature.csv) (human to mouse converted) 
2. regenerated in /Users/gs/G/data/MK_myc_2022/myc_mouse/scripts/01_load_data.R with:
    - added custom PRO and ANTI-apoptotic sets
    - mouse MSigDB Hallmark pathways.
3. Saved as saveRDS(gene_sets_list, here("results", "gene_sets_list.rds")
This list was used throughout pathway analysis in the myc_mouse/main branch scripts 04_fgsea... 07_heatmap, 08_pathway_summary.
Then:
4. Cell death genes:
    - /Users/gs/G/data/MK_myc_2022/myc_mouse/data/cell_death_genes_consolidated.rds this is a compilation of cell death related genes from GO, KEGG, Reactome, MSIgDB - categorised by pathway (Apoptosis, CICD), pro vs anti, mitochondrial and ''core''. Contin both mouse and human ENSEMBL and gene symbols. Used MyGene.info API with HomoloGene database for ortholog lookup (93% ). Was used on hypothesis testing whether pro-death effect is reduced at 12W. /Users/gs/G/data/MK_myc_2022/myc_mouse/scripts/archive_main_pipeline/09_cell_death_pathway_summary.R

**Myc-mouse/new_analysis**

5. Overall pathway analysis uses the same /Users/gs/G/data/MK_myc_2022/myc_mouse/results/gene_sets_list.rds as the main branch.
6. MitoPPS analysis (/Users/gs/G/data/MK_myc_2022/myc_mouse/scripts/08_mitoPPS_analysis.R) uses full mitocarta3.0 pathways, modified: mtDNA encoded OXPHOS subunits are separated from OXPHOS (mtOXPHOS and nuOXPHOS, Pro and Anti-apoptotic separation in APOPTOTIC mitocarta pathway). fGSEA was also used in this analysis branch (/Users/gs/G/data/MK_myc_2022/myc_mouse/scripts/09_mitoPPS_vs_fgsea_comparison.R). THen also used in interaction analysis.
7. Cell death pathay analysis uses different gene sets. 15 gene sets from [Tang et al](https://doi.org/10.1016/j.csbj.2024.08.012)  (/Users/gs/G/data/MK_myc_2022/myc_mouse/data/cell-death) 

based on this a [[prompt for the  mammary_geneset_library to amend claude code plan]]