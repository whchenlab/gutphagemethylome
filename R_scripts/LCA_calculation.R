calculate_LCA<-function(host2virus){
  host2virus<-data.frame(host2virus)
  colnames(host2virus)<-c("Virus","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  host2virus$Species<-gsub(" NA","",host2virus$Species)
  viruscount<-host2virus$Virus %>% table%>% data.frame()
  one<-subset(viruscount,Freq==1)[,1]
  more<-subset(viruscount,Freq>1)[,1]
  sub1<-subset(host2virus,Virus %in% one & Species!="s__")[,c("Species","Virus")]
  sub2list<-subset(host2virus,Virus %in% one & Species=="s__")[,"Virus"]
  sub2<-subset(host2virus,Virus %in% sub2list)[,c("Species","Virus")]
  # sub2$Host<-gsub("s__ NA","s__",sub2$Host)
  names(sub2)<-names(sub1)
  subone<-rbind(sub1,sub2)
  if(!length(subone$Species)==0){
    subone$"level"<-"Species"
    names(subone)<-c("LCA","Virus","level")
  }
  submore<- data.frame(LCA="LCA",Virus="Virus",level="level")
  for(i in more){
    link<-subset(host2virus,Virus==i)
    for(k in 2:8){
      num<-unique(link[,k]) %>% length
      if(num>1 & k==2){
        level="Kingdom"
        name<-""
        break
      }
      if(num>1){
        level=switch(k-2,"Kingdom","Phylum","Class","Order","Family","Genus","Species")
        name<-unique(link[,k-1])
        if(unique(link[,k])=="s__"){
          level=switch(k-2,"Kingdom","Phylum","Class","Order","Family","Genus","Species")
          name=unique(link[,k-1])
        }
        break
      }
      #print(i)
      if(num==1 & k==8){
        level="Species"
        name<-unique(link[,k])
        if(unique(link[,k])=="s__"){
          level="Genus"
          name=unique(link[,k-1])
        }
        break
      }
    }
    lca<-c(name,i,level)
    submore<-rbind(submore,lca)
  }
  #names(submore)<-c("LCA","Virus","level")
  submore<-submore[-1,]
  all<-rbind(subone,submore)
  return(all)
}
