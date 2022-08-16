## Shell scripts
> Rawdatas deposited in [https://figshare.com/articles/dataset/Rawdata/20464083](https://figshare.com/articles/dataset/Rawdata/20464203)

> CHGV.hq.faa -- predicted CHGV proteins
> 
> CHGV.hq.MTase.faa -- predicted CHGV MTases
> 
> HGV.hq.genome.fa -- CHGV HQ viral genomes
> 
> UHGGv2.MTase.faa -- predicted UHGG MTases
> 
> Modified_gff.rar -- basemod gff files generated from SMRTlink
### 1. Batch_upload_to_SMRTlink.sh
> upload xml files to SMRTlink websever with import-dataset API
### 2. Dealing_with_modified_gffs.sh
> Integrate SMRTlink generated basemods gff files. 
> 
> Motif discovery and locating
### 3. get_MTases.sh
> Search for MTases in selected proteins with conserved domain.
### 4. MTase.blast.mcl.sh
> All-against-all blast of CHGV MTases and UHGGv2-UHGP95 MTases.
> 
> MCL clusteing using their similarity.

## Rscripts
> R datas deposited in https://figshare.com/articles/dataset/Rdatas/20464083
### 1. Dealing_with_motifs.R
> Functions to remove redundancy and verify that whether the selected motifs exist in REbase.
### 2. Figures
> Rscripts to plot the main figures.
### 3. load_functions.R
```r
calculate_LCA() # custom function to calculate the last common anscestor of selected bacterias
get_UP_degenerate_bases()
get_degenerate_bases() # functions to find the range of selected motis
```
