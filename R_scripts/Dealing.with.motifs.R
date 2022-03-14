library(Biostrings)
library(dplyr)
motif.list<-motif.filtered[,c(1,5,8)] %>%unique

get_UP_degenerate_bases<-function(s){
  ss<-gsub("D","[DN]",s)
  ss<-gsub("V","[VN]",ss)
  ss<-gsub("B","[BN]",ss)
  ss<-gsub("H","[HN]",ss)
  ss<-gsub("W","[WHDN]",ss)
  ss<-gsub("S","[SBVN]",ss)
  ss<-gsub("K","[KBDN]",ss)
  ss<-gsub("M","[MHVN]",ss)
  ss<-gsub("Y","[YHBN]",ss)
  ss<-gsub("R","[RVDN]",ss)
  ss<-gsub("C","[CYMSHBVN]",ss)
  ss<-gsub("G","[GRKSBVDN]",ss)
  ss<-gsub("T","[TYKWHBDN]",ss)
  ss<-gsub("A","[ARMWHVDN]",ss)
  return(ss)
}
get_degenerate_bases<-function(s){
  ss<-gsub("R","[AGR]",s)
  ss<-gsub("Y","[CTY]",ss)
  ss<-gsub("M","[ACM]",ss)
  ss<-gsub("K","[GTK]",ss)
  ss<-gsub("S","[GCS]",ss)
  ss<-gsub("W","[ATW]",ss)
  ss<-gsub("H","[ATCMWYH]",ss)
  ss<-gsub("B","[GTCKYSB]",ss)
  ss<-gsub("V","[GACRMSV]",ss)
  ss<-gsub("D","[GATRKWD]",ss)
  ss<-gsub("N","[ATCGRYMKSWHBVDN]",ss)
  return(ss)
}
motif.list$exist<-NA
for(s in (MTase.blast.filtered$RecSeq %>% unique %>% na.omit) ){
  #motif.list[grep(get_degenerate_bases(s),motif.list$motif_id)]
  if(length(grep(get_degenerate_bases(s),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(s),motif.list$motif_id)){
    start<-regexpr(get_degenerate_bases(s), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(s), motif.list$motif_id[i])[1]+nchar(s)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$exist[i]<-s
    }
  }
  ss<-stringi::stri_reverse(s)
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
    start<-regexpr(ss, motif.list$motif_id[i])[1]
    stop<-regexpr(ss, motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$exist[i]<-s
    }
  }
  ss<-Biostrings::complement(DNAString(s)) %>% as.character()
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
    start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$exist[i]<-s
    }
  }
  ss<-Biostrings::reverseComplement(DNAString(s)) %>% as.character()
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
    start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$exist[i]<-s
    }
  }
}
motif.list$partial<-NA
for(s in (MTase.blast.filtered$RecSeq %>% unique %>% na.omit) ){
  #motif.list[grep(get_degenerate_bases(s),motif.list$motif_id)]
  if(length(grep(get_UP_degenerate_bases(s),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_UP_degenerate_bases(s),motif.list$motif_id)){
    start<-regexpr(get_UP_degenerate_bases(s), motif.list$motif_id[i])[1]
    stop<-regexpr(get_UP_degenerate_bases(s), motif.list$motif_id[i])[1]+nchar(s)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$partial[i]<-s
    }
  }
  ss<-stringi::stri_reverse(s)
  if(length(grep(get_UP_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_UP_degenerate_bases(ss),motif.list$motif_id)){
    start<-regexpr(get_UP_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_UP_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$partial[i]<-s
    }
  }
  ss<-Biostrings::complement(DNAString(s)) %>% as.character()
  if(length(grep(get_UP_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_UP_degenerate_bases(ss),motif.list$motif_id)){
    start<-regexpr(get_UP_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_UP_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$partial[i]<-s
    }
  }
  ss<-Biostrings::reverseComplement(DNAString(s)) %>% as.character()
  if(length(grep(get_UP_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_UP_degenerate_bases(ss),motif.list$motif_id)){
    start<-regexpr(get_UP_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_UP_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$partial[i]<-s
    }
  }
}

motif.list$partial[which(!is.na(motif.list$partial)&!is.na(motif.list$exist))]<-NA

#motif dereplication----
motif.list$same<-NA
for(s in (motif.list$motif_id %>% unique %>% na.omit) ){
  #motif.list[grep(get_degenerate_bases(s),motif.list$motif_id)]
  if(length(grep(get_degenerate_bases(s),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(s),motif.list$motif_id)){
    if((!is.na(motif.list$same[i])&&nchar(s)>nchar(motif.list$same[i]))||s==motif.list$motif_id[i]){
      break
    }
    start<-regexpr(get_degenerate_bases(s), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(s), motif.list$motif_id[i])[1]+nchar(s)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$same[i]<-s
    }
  }
  ss<-stringi::stri_reverse(s)
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
    if((!is.na(motif.list$same[i])&&nchar(s)>nchar(motif.list$same[i]))||s==motif.list$motif_id[i]){
      break
    }
    start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$same[i]<-s
    }
  }
  ss<-Biostrings::complement(DNAString(s)) %>% as.character()
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
    if((!is.na(motif.list$same[i])&&nchar(s)>nchar(motif.list$same[i]))||s==motif.list$motif_id[i]){
      break
    }
    start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$same[i]<-s
    }
  }
  ss<-Biostrings::reverseComplement(DNAString(s)) %>% as.character()
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))==0){
    next
  }
  for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
    if((!is.na(motif.list$same[i])&&nchar(s)>nchar(motif.list$same[i]))||s==motif.list$motif_id[i]){
      break
    }
    start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
    stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
    if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
      motif.list$same[i]<-s
    }
  }
}

motif.list$same[which(is.na(motif.list$same))]<-motif.list$motif_id[which(is.na(motif.list$same))]
