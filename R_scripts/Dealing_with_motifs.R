motif.list<-motif.filtered[,c(1,2,8)] %>%unique
motif.list$exist<-NA
for(s in (RE.info[grep("methyltransferase",RE.info$Type),"RecSeq"] %>% unique ) ){
  if(length(grep(s,motif.list$exist))!=0){
    next
  }
  if(is.na(s)){
    next
  }
  if(length(grep(0,s))!=0){
    next
  }
  if(length(grep(get_degenerate_bases(s),motif.list$motif_id))!=0){
    for(i in grep(get_degenerate_bases(s),motif.list$motif_id)){
      start<-regexpr(get_degenerate_bases(s), motif.list$motif_id[i])[1]
      stop<-regexpr(get_degenerate_bases(s), motif.list$motif_id[i])[1]+nchar(s)-1
      if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
        motif.list$exist[i]<-s
      }
    }
  }
  
  ss<-stringi::stri_reverse(s)
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))!=0){
    for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
      start<-regexpr(ss, motif.list$motif_id[i])[1]
      stop<-regexpr(ss, motif.list$motif_id[i])[1]+nchar(ss)-1
      if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
        motif.list$exist[i]<-s
      }
    }
  }
  
  ss<-Biostrings::complement(DNAString(s)) %>% as.character()
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))!=0){
    for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
      start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
      stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
      if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
        motif.list$exist[i]<-s
      }
    }
  }
  
  ss<-Biostrings::reverseComplement(DNAString(s)) %>% as.character()
  if(length(grep(get_degenerate_bases(ss),motif.list$motif_id))!=0){
    for(i in grep(get_degenerate_bases(ss),motif.list$motif_id)){
      start<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]
      stop<-regexpr(get_degenerate_bases(ss), motif.list$motif_id[i])[1]+nchar(ss)-1
      if(motif.list$position[i]>=start && motif.list$position[i]<=stop){
        motif.list$exist[i]<-s
      }
    }
  }
  
}
