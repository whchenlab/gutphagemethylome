## Data

> Related data files are deposited to [https://figshare.com/articles/dataset/Rawdata/20464083](https://figshare.com/articles/dataset/Rawdata/20464203)

### Sequence files 

> HGV.hq.genome.fa -- CHGV HQ viral genomes
> 
> CHGV.hq.faa -- predicted CHGV proteins
> 
> CHGV.hq.MTase.faa -- predicted CHGV MTases
> 
> UHGGv2.MTase.faa -- predicted UHGG MTases
> 
> pt2genome.csv -- protein id to genome name
> 
### Processed methylome data
> Modified_gff.rar -- basemod gff files generated from SMRTlink for the CHGV HQ viral genomes; one GFF file for each sample.

### Rdata files

> various RData files required for the Rscripts below.

## Shell scripts
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
