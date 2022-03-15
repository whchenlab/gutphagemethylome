#load packages----
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(forcats)
library(RColorBrewer)
library(gg.gap)
library(patchwork)
library(reshape2)
library(pROC)

#load data----
load("Data/modified.anno.cds.RData")
load("Data/genome.stat.RData")
load("Data/lifestyle.RData")
load("Data/Fun.RData")
load("Data/genome.anno.info.RData")
load("Data/gene.stat.RData")
load("Data/motif.filtered.RData")
load("Data/genome.info.Rdata")
load("Data/genome.gene.count.Rdata")
load("Data/GPD.meta.RData")
load("Data/ases.blast.RData")
load("Data/mvp.blast.RData")

#Figure1C----
t<-unique(modified.anno.cds.done[,c(1,3,4,10)])
modified.anno.cds.done.stat<-aggregate(t$motif_site,
                                       by=list(t$sequence_name,
                                               t$InCDS,
                                               t$Modified_Type),
                                       FUN=length)
names(modified.anno.cds.done.stat)<-c("sequence_name","InCDS","Modified_Type","Modified_Bases")
modified.anno.cds.done.stat$Group<-"In CDS"
modified.anno.cds.done.stat$Group[which(modified.anno.cds.done.stat$InCDS=="N")]<-"Non-coding region"
modified.anno.cds.done.stat<-merge(modified.anno.cds.done.stat,genome.stat,
                                   by.x="sequence_name",
                                   by.y="Genome")
modified.anno.cds.done.stat$Modified_Bases_pct<-NA
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m6A"&modified.anno.cds.done.stat$InCDS=="Y")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m6A"& modified.anno.cds.done.stat$InCDS=="Y")]/modified.anno.cds.done.stat$AT.cds[which(modified.anno.cds.done.stat$Modified_Type=="m6A"& modified.anno.cds.done.stat$InCDS=="Y")]
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m6A"&modified.anno.cds.done.stat$InCDS=="N")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m6A"& modified.anno.cds.done.stat$InCDS=="N")]/modified.anno.cds.done.stat$nonCDS.AT[which(modified.anno.cds.done.stat$Modified_Type=="m6A"&modified.anno.cds.done.stat$InCDS=="N")]
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="Y")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m4C"& modified.anno.cds.done.stat$InCDS=="Y")]/modified.anno.cds.done.stat$GC.cds[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="Y")]
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="N")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m4C"& modified.anno.cds.done.stat$InCDS=="N")]/modified.anno.cds.done.stat$nonCDS.GC[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="N")]

modified.anno.cds.done.stat %>% 
  ggplot(aes(x=Group,y=Modified_Bases_pct))+
  geom_boxplot(aes(color=Group),width=0.5)+
  scale_y_log10()+
  theme_bw()+
  scale_color_futurama()+
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat$Group)),step_increase = 0.1,
              map_signif_level = T,test = t.test)+
  labs(y="%modified bases")+
  facet_wrap(~Modified_Type)+
  theme(strip.background.x = element_rect(fill = "white")) +
  guides(color="none")


t<-unique(modified.anno.cds.done[,c(1,3,4,10)])
modified.anno.cds.done.stat<-aggregate(t$motif_site,
                                       by=list(t$sequence_name,
                                               t$InCDS,
                                               t$Modified_Type),
                                       FUN=length)
names(modified.anno.cds.done.stat)<-c("sequence_name","InCDS","Modified_Type","Modified_Bases")
modified.anno.cds.done.stat$Group<-"In CDS"
modified.anno.cds.done.stat$Group[which(modified.anno.cds.done.stat$InCDS=="N")]<-"Non-coding region"
modified.anno.cds.done.stat<-merge(modified.anno.cds.done.stat,genome.stat,
                                   by.x="sequence_name",
                                   by.y="Genome")
modified.anno.cds.done.stat$Modified_Bases_pct<-NA
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m6A"&modified.anno.cds.done.stat$InCDS=="Y")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m6A"& modified.anno.cds.done.stat$InCDS=="Y")]/modified.anno.cds.done.stat$AT.cds[which(modified.anno.cds.done.stat$Modified_Type=="m6A"& modified.anno.cds.done.stat$InCDS=="Y")]
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m6A"&modified.anno.cds.done.stat$InCDS=="N")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m6A"& modified.anno.cds.done.stat$InCDS=="N")]/modified.anno.cds.done.stat$nonCDS.AT[which(modified.anno.cds.done.stat$Modified_Type=="m6A"&modified.anno.cds.done.stat$InCDS=="N")]
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="Y")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m4C"& modified.anno.cds.done.stat$InCDS=="Y")]/modified.anno.cds.done.stat$GC.cds[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="Y")]
modified.anno.cds.done.stat$Modified_Bases_pct[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="N")] <- modified.anno.cds.done.stat$Modified_Bases[which(modified.anno.cds.done.stat$Modified_Type=="m4C"& modified.anno.cds.done.stat$InCDS=="N")]/modified.anno.cds.done.stat$nonCDS.GC[which(modified.anno.cds.done.stat$Modified_Type=="m4C"&modified.anno.cds.done.stat$InCDS=="N")]

modified.anno.cds.done.stat<-merge(modified.anno.cds.done.stat,lifestyle[,c(1,4)],by.x="sequence_name",by.y="X1")
modified.anno.cds.done.stat %>% 
  ggplot(aes(x=Group,y=Modified_Bases_pct))+
  geom_boxplot(aes(color=Group),width=0.5)+
  scale_y_log10(labels = scales::percent)+
  theme_bw()+
  scale_color_futurama()+
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat$Group)),step_increase = 0.1,
              map_signif_level = T,test = t.test)+
  labs(y="%modified bases")+
  facet_grid(Modified_Type~X4)+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white")) +
  guides(color="none")

#Figure1D----
names(anno)<-c("sequence_name","pt_id","Func","Func_group")
anno<-subset(anno,sequence_name %in% modified.anno.cds.done.stat$sequence_name)
anno<-unique(anno)
names(genome.anno.info)<-c("sequence_name","Type","start","stop","strand","pt_id","product")
anno<-merge(anno,genome.anno.info,all=T,by=c("sequence_name","pt_id"))
anno$Func_group[which(is.na(anno$Func_group))]<-"Hypothetical protein"
anno.stat<-aggregate(unique(anno[,c(1,2)])$pt_id,by=list(unique(anno[,c(1,2)])$sequence_name),FUN=length)
names(anno.stat)<-c("sequence_name","pt_count")
t1<-modified.anno.cds.done[,c(1,2,4)] %>% unique
t1<-subset(t1,!is.na(pt_id))
modified.anno.cds.done.stat.stat<-aggregate(t1$pt_id,
                                            by=list(t1$sequence_name,
                                                    t1$Modified_Type),
                                            FUN=length)
names(modified.anno.cds.done.stat.stat)<-c("sequence_name","Modified_Type","Modified_pts")
modified.anno.cds.done.stat.stat<-merge(modified.anno.cds.done.stat.stat,anno.stat,by=c("sequence_name"))
modified.anno.cds.done.stat.stat<-merge(modified.anno.cds.done.stat.stat,lifestyle[,c(1,4)],by.x="sequence_name",by.y="X1")

modified.anno.cds.done.stat.stat %>% 
  ggplot(aes(x=Modified_Type,y=Modified_pts/pt_count))+
  geom_boxplot(aes(fill=Modified_Type),width=0.5,outlier.size = 0.1)+
  scale_fill_manual(values=c("#bce7cb","#edcbc2"))+
  theme_bw()+
  guides(fill="none")+
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat.stat$Modified_Type)),step_increase = 0.1,
              map_signif_level = T,test = wilcox.test)+
  
  labs(y="% modified",x="Methylation type")+
  theme(strip.background.x = element_rect(fill = "white"),axis.text.x = element_blank()) 

modified.anno.cds.done.stat.stat %>% 
  ggplot(aes(x=Modified_Type,y=Modified_pts/pt_count))+
  geom_boxplot(aes(fill=Modified_Type),width=0.5,outlier.size = 0.1)+
  scale_fill_manual(values=c("#bce7cb","#edcbc2"))+
  theme_bw()+
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat.stat$Modified_Type)),step_increase = 0.1,
              map_signif_level = T,test = wilcox.test)+
  facet_grid(~X4)+
  labs(y="% modified",x="Methylation type")+
  theme(strip.background.x = element_rect(fill = "white"),axis.text.x = element_blank()) 

#Figure1E----
load("Data/Fun.RData")

names(anno)<-c("sequence_name","pt_id","Func","Func_group")

anno<-subset(anno,sequence_name %in% modified.anno.cds.done.stat$sequence_name)
anno<-unique(anno)
names(genome.anno.info)<-c("sequence_name","Type","start","stop","strand","pt_id","product")
anno<-merge(anno,genome.anno.info,all=T,by=c("sequence_name","pt_id"))
anno$Func_group[which(is.na(anno$Func_group))]<-"Hypothetical protein"
anno.stat<-aggregate(unique(anno[,c(1,2,4)])$pt_id,by=list(unique(anno[,c(1,2,4)])$sequence_name,unique(anno[,c(1,2,4)])$Func_group),FUN=length)
names(anno.stat)<-c("sequence_name","Function_group","pt_count")

modified.anno.cds.done.stat.t<-unique(modified.anno.cds.done[,c("Modified_Type","sequence_name","Func_group","pt_id")])
modified.anno.cds.done.stat.stat<-aggregate(modified.anno.cds.done.stat.t$Modified_Type,
                                            by=list(modified.anno.cds.done.stat.t$sequence_name,
                                                    modified.anno.cds.done.stat.t$Modified_Type,
                                                    modified.anno.cds.done.stat.t$Func_group),
                                            FUN=length)

names(modified.anno.cds.done.stat.stat)<-c("sequence_name","Modified_Type","Function_group","Modified_pts")
modified.anno.cds.done.stat.stat<-merge(modified.anno.cds.done.stat.stat,anno.stat,by=c("sequence_name","Function_group"))
modified.anno.cds.done.stat.stat$"ng"<-"Known Viral function"
modified.anno.cds.done.stat.stat$"ng"[which(modified.anno.cds.done.stat.stat$Function_group=="tRNA")]<-"tRNA"
modified.anno.cds.done.stat.stat$"ng"[which(modified.anno.cds.done.stat.stat$Function_group=="Unsorted")]<-"Unsorted"
modified.anno.cds.done.stat.stat$"ng"[which(modified.anno.cds.done.stat.stat$Function_group=="Hypothetical protein")]<-"Hypothetical protein"

modified.anno.cds.done.stat.stat %>% 
  mutate(ng=fct_reorder(ng,-(Modified_pts/pt_count),.fun = mean)) %>%
  ggplot(aes(x=Modified_Type,y=Modified_pts/pt_count))+
  geom_boxplot(aes(fill=Modified_Type),width=0.5,outlier.size = 0.1)+
  scale_fill_manual(values=c("#EF7548","#77A577"))+
  theme_bw()+
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat.stat$Modified_Type)),step_increase = 0.1,
              map_signif_level = T,test = wilcox.test)+
  facet_wrap(~ng,nrow=1)+
  scale_y_continuous(labels=scales::percent)+

  labs(y="% modified CDS")+
  theme(strip.background.x = element_rect(fill = "white"),axis.text.x = element_blank()) 

#Figure1F----
modified.anno.cds.done.stat.t<-unique(modified.anno.cds.done[,c("motif_site","sequence_name","pt_id","Modified_Type","Func_group")])
modified.anno.cds.done.stat.stat<-aggregate(modified.anno.cds.done.stat.t$motif_site,
                                            by=list(modified.anno.cds.done.stat.t$sequence_name,
                                                    modified.anno.cds.done.stat.t$pt_id,
                                                    modified.anno.cds.done.stat.t$Modified_Type,
                                                    modified.anno.cds.done.stat.t$Func_group),
                                            FUN=length)
names(modified.anno.cds.done.stat.stat)<-c("sequence_name","protein id","Modified_Type","Function_group","Modified_bases")
modified.anno.cds.done.stat.stat<-merge(modified.anno.cds.done.stat.stat,lifestyle[,c(1,4)],by.x="sequence_name",by.y="X1")
gene.stat.1<-gene.stat[,c(1,2,8,9)]
colnames(gene.stat.1)<-c("sequence_name","protein id","GC","AT")
modified.anno.cds.done.stat.stat<-merge(modified.anno.cds.done.stat.stat,gene.stat.1,by=c("sequence_name","protein id"))
modified.anno.cds.done.stat.stat$Total<-NA
modified.anno.cds.done.stat.stat$Total[which(modified.anno.cds.done.stat.stat$Modified_Type=="m6A")] <- modified.anno.cds.done.stat.stat$AT[which(modified.anno.cds.done.stat.stat$Modified_Type=="m6A")]
modified.anno.cds.done.stat.stat$Total[which(modified.anno.cds.done.stat.stat$Modified_Type=="m4C")] <- modified.anno.cds.done.stat.stat$GC[which(modified.anno.cds.done.stat.stat$Modified_Type=="m4C")]

modified.anno.cds.done.stat.stat %>% 
  subset(sequence_name %in% (modified.anno.cds.done.stat.t$sequence_name %>% unique)) %>% 
  mutate(Function_group=fct_reorder(Function_group,-Modified_bases)) %>%
  ggplot(aes(x=Modified_Type,y=Modified_bases/Total))+
  geom_boxplot(aes(fill=Modified_Type),width=0.5,outlier.size = 0.1)+
  scale_y_log10()+
  theme_bw()+
  scale_fill_manual(values=c("#bce7cb","#edcbc2"))+
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat.stat$Modified_Type)),step_increase = 0.1,
              map_signif_level = T,test = wilcox.test)+
  guides(fill="none")+
  labs(y="%modified bases per CDS ",x="Methylation type")+
  theme(strip.background.x = element_rect(fill = "white"),axis.text.x = element_blank()) 
modified.anno.cds.done.stat.stat %>% 
  subset(sequence_name %in% (modified.anno.cds.done.stat.t$sequence_name %>% unique)) %>% 
  mutate(Function_group=fct_reorder(Function_group,-Modified_bases)) %>%
  ggplot(aes(x=Modified_Type,y=Modified_bases/Total))+
  geom_boxplot(aes(fill=Modified_Type),width=0.5,outlier.size = 0.1)+
  scale_y_log10()+
  theme_bw()+
  scale_fill_manual(values=c("#bce7cb","#edcbc2"))+
  
  geom_signif(comparisons = list(unique(modified.anno.cds.done.stat.stat$Modified_Type)),step_increase = 0.1,
              map_signif_level = T,test = wilcox.test)+
  facet_wrap(~Function_group,nrow=2)+
  # facet_grid(~X4)+
  labs(y="%modified bases per CDS ",x="Methylation type")+
  theme(strip.background.x = element_rect(fill = "white"),axis.text.x = element_blank()) 


#Figure1H----
modified.anno.cds.done.1<-modified.anno.cds.done[,c(1,2,3,4,10)]
colnames(modified.anno.cds.done.1)[4]<-"Type"
motif.filtered.pt<-merge(motif.filtered[,c(3,4,5,9)],modified.anno.cds.done.1,by=c("sequence_name","Type","motif_site"))
motif.filtered.pt.1<-motif.filtered.pt[,c(1,2,4,6)] %>% unique
motif.filtered.pt.2<-motif.filtered.pt[,c(2,4,6)] %>% unique
motif.filtered.pt.agg<-aggregate(motif.filtered.pt.1$same,
                                 by=list(motif.filtered.pt.1$sequence_name,motif.filtered.pt.1$Type,motif.filtered.pt.1$InCDS),
                                 FUN=length)
colnames(motif.filtered.pt.agg)<-c("sequence_name","Type","InCDS","count")
motif.filtered.pt.agg<-merge(motif.filtered.pt.agg,lifestyle[,c(1,4)],by.x="sequence_name",by.y="X1")
motif.filtered.pt.agg.1<-aggregate(motif.filtered.pt.2$same,
                                   by=list(motif.filtered.pt.2$Type,motif.filtered.pt.2$InCDS),
                                   FUN=length)
m4C=(subset(motif.filtered.pt,Type=="m4C")$same %>% na.omit() %>% unique %>% length)
m6A=(subset(motif.filtered.pt,Type=="m6A")$same %>% na.omit() %>% unique %>% length)
m4C.cds<-(motif.filtered.pt.agg.1 %>% subset(Group.1=="m4C"&Group.2=="Y"))$x
m4C.ncds<-(motif.filtered.pt.agg.1 %>% subset(Group.1=="m4C"&Group.2=="N"))$x
m6A.cds<-(motif.filtered.pt.agg.1 %>% subset(Group.1=="m6A"&Group.2=="Y"))$x
m6A.ncds<-(motif.filtered.pt.agg.1 %>% subset(Group.1=="m6A"&Group.2=="N"))$x
m4C.both<-m4C.cds+m4C.ncds-m4C
m6A.both<-m6A.cds+m6A.ncds-m6A
p0<-data.frame(Type=c(rep("m4C",3),rep("m6A",3)),
               Group=rep(c("CDS Region","Both","Non-coding region"),2),
               value=c(m4C.cds-m4C.both,m4C.both,m4C.ncds-m4C.both,m6A.cds-m6A.both,m6A.both,m6A.ncds-m6A.both)) %>% 
  mutate(Group=factor(Group,levels=c("CDS Region","Both","Non-coding region"))) %>% 
  ggplot(aes(x=Type,y=value))+
  geom_bar(aes(fill=Group),stat = "identity") +
  scale_fill_manual(values=c("#da9959","#bc6c35","#f7deba"))+
  theme_bw()+
  guides(fill="none")+
  labs(x="Methylation Type",y="Nr. of Motif Type")

p00<-gg.gap(plot = p0,
            segments = c(10, 300),
            tick_width = 50,
            rel_heights = c(0.05, 0, 0.1),
            ylim = c(0, 700))
p1<-motif.filtered.pt.agg%>% 
  # mutate(Gene_type=factor(Gene_type,levels = c("MTase","REase"))) %>% 
  # mutate(count=factor(count,levels = c("0","1","2","3","4",">=5"))) %>% 
  ggplot(aes(x=InCDS,y=count,fill=InCDS))+
  geom_boxplot(width=0.5)+
  scale_fill_brewer(palette = "GnBu")+
  theme_bw()+
  labs(x="")+
  geom_signif(comparisons = list(unique(motif.filtered.pt.agg$InCDS)),step_increase = 0.05,
              map_signif_level = T,test = t.test)+
  facet_grid(~Type)+
  theme(axis.text.x = element_blank())+
  scale_y_log10()+
  theme(strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_rect(fill = "white"),
        axis.text.x = element_blank()) +
  labs(y="Nr. of Motif Type",color="InCDS")
p2<-motif.filtered.pt.agg%>% 
  # mutate(Gene_type=factor(Gene_type,levels = c("MTase","REase"))) %>% 
  # mutate(count=factor(count,levels = c("0","1","2","3","4",">=5"))) %>% 
  ggplot(aes(x=InCDS,y=count,fill=InCDS))+
  geom_boxplot(width=0.5)+
  geom_signif(comparisons = list(unique(motif.filtered.pt.agg$InCDS)),step_increase = 0.05,
              map_signif_level = T,test = t.test)+
  scale_fill_brewer(palette = "GnBu")+
  # scale_color_aaas()+
  theme_bw()+
  labs(x="")+
  facet_grid(Type~X4)+
  theme(strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_rect(fill = "white"),
        axis.text.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  scale_y_log10()+
  labs(y="Nr. of Motif Type",color="InCDS")
p1.1<-p1+p2+plot_layout(width=c(1,2),guides = 'collect')
p<-p00+p1.1+plot_layout(width=c(1,4),guides = 'collect')
p


#Figure2B----
t.stat<-aggregate(t$motif_site,by=list(t$sequence_name,t$Modified_Type),FUN=length)
names(t.stat)<-c("sequence_name","Modified_Type","Modified_Bases")
t.stat<-merge(t.stat,genome.stat,
              by.x="sequence_name",
              by.y="Genome")
t.stat$Modified_Bases_pct<-0
t.stat$Modified_Bases_pct[which(t.stat$Modified_Type=="m6A")] <- t.stat$Modified_Bases[which(t.stat$Modified_Type=="m6A")]/t.stat$AT[which(t.stat$Modified_Type=="m6A")]
t.stat$Modified_Bases_pct[which(t.stat$Modified_Type=="m4C")] <- t.stat$Modified_Bases[which(t.stat$Modified_Type=="m4C")]/t.stat$GC[which(t.stat$Modified_Type=="m4C")]

genome.info.count<-merge(genome.info,genome.gene.count,by.x="names",by.y="Genome",all.x=T)
genome.info.count<-merge(genome.info.count,t.stat,by.x="names",by.y="sequence_name")

genome.info.count$count[which(is.na(genome.info.count$count))]<-0
r<-genome.info.count[which(is.na(genome.info.count$Gene_type)),]
r$Gene_type<-"REase"
genome.info.count$Gene_type[which(is.na(genome.info.count$Gene_type))]<-"MTase"
genome.info.count<-rbind(genome.info.count,r)
genome.info.count$count[which(genome.info.count$count>=5)]<-">=5"
genome.info.count<-genome.info.count[,c("names","Gene_type","count","Modified_Type","Modified_Bases_pct")]
genome.info.count<-dcast(genome.info.count,names+Modified_Type+Modified_Bases_pct~Gene_type,value.var = "count")
genome.info.count[is.na(genome.info.count)]<-0
genome.info.count<-melt(genome.info.count,
                        id=c(c("names","Modified_Type","Modified_Bases_pct")),
                        variable.name = "Gene_type",
                        value.name = "count")
save(genome.info.count,file="../genome.info.count.Rdata")
genome.info.count %>% 
  mutate(Gene_type=factor(Gene_type,levels = c("MTase","REase"))) %>% 
  mutate(count=factor(count,levels = c("0","1","2","3","4",">=5"))) %>% 
  ggplot(aes(x=count,y=Modified_Bases_pct,color=Modified_Type))+
  geom_boxplot(width=0.6)+
  scale_color_futurama()+
  theme_bw()+
  scale_y_log10(labels = scales::percent)+
  facet_grid(~Gene_type)+
  # theme(scales=scales::percent)+
  labs(x="Gene count",y="% bases modified",color="Modified_Type")+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))

# motif.filtered<-motif.filtered[,c(1,2,3,4,5,6,7,8)] %>% unique
# motif.filtered<- merge(motif.filtered,motif.list[,c(1,2,3,6)],by=c("motif_id","position","Type")) 
motif.filtered.count<-aggregate(unique(motif.filtered[,c(3,4,9)])$same,
                                by=list(unique(motif.filtered[,c(3,4,9)])$sequence_name,
                                        unique(motif.filtered[,c(3,4,9)])$Type),FUN=length)
names(motif.filtered.count)<-c("Genome","Motif_Type","motif_count")


genome.info.count<-merge(genome.info,genome.gene.count,by.x="names",by.y="Genome",all.x=T)
genome.info.count$count[which(is.na(genome.info.count$count))]<-0

r<-genome.info.count[which(is.na(genome.info.count$Gene_type)),]
r$Gene_type<-"REase"
genome.info.count$Gene_type[which(is.na(genome.info.count$Gene_type))]<-"MTase"
genome.info.count<-rbind(genome.info.count,r)
genome.info.count$count[which(genome.info.count$count>=5)]<-">=5"
genome.info.count<-genome.info.count[,c("names","Gene_type","count")]
genome.info.count<-dcast(genome.info.count,names~Gene_type,value.var = "count")
genome.info.count[is.na(genome.info.count)]<-0
genome.info.count<-melt(genome.info.count,
                        id=c(c("names")),
                        variable.name = "Gene_type",
                        value.name = "count")

genome.info.count<-merge(genome.info.count,motif.filtered.count,by.x="names",by.y="Genome",all.x=T)
genome.info.count$count[which(genome.info.count$count>=5)]<-">=5"
genome.info.count %>% 
  mutate(Gene_type=factor(Gene_type,levels = c("MTase","REase"))) %>% 
  mutate(count=factor(count,levels = c("0","1","2","3","4",">=5"))) %>% 
  ggplot(aes(x=count,y=motif_count,fill=Motif_Type))+
  geom_boxplot(width=0.6)+
  scale_fill_manual(values=c("#ef754b","#77a577"))+
  theme_bw()+
  facet_grid(~Gene_type)+
  scale_y_log10()+
  labs(x="Gene count",y="Nr. of Motif",fill="Motif Type")+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))
#Figure2C----
chi.test.data<-aggregate(genome.info$MTase,by=list(genome.info$Lifestyle,genome.info$MTase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chisq.test(chi.test.data[,c(2,3)])
chisq.test(chi.test.data[,c(3,4)])
chisq.test(chi.test.data[,c(4,5)])
mt<-genome.info %>% 
  mutate(MTase=factor(MTase,levels = c("Y","N"))) %>% 
  ggplot()+
  geom_bar(aes(x=Lifestyle,y=..count..,fill=MTase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+
  labs(y="Percentage",title = "MTase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))

chi.test.data<-aggregate(genome.info$REase,by=list(genome.info$Lifestyle,genome.info$REase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chisq.test(chi.test.data[,c(2,5)])
chisq.test(chi.test.data[,c(2,3)])
chisq.test(chi.test.data[,c(3,4)])
chisq.test(chi.test.data[,c(4,5)])
re<-genome.info %>% 
  mutate(REase=factor(REase,levels = c("Y","N"))) %>% 
  ggplot()+
  geom_bar(aes(x=Lifestyle,y=..count..,fill=REase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+
labs(y="",title = "REase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))

chi.test.data<-aggregate(genome.info$MTase,by=list(genome.info$Type,genome.info$MTase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chi.test.data[is.na(chi.test.data)]<-0
chisq.test(chi.test.data[,c("High-quality","crAssphage")])
chisq.test(chi.test.data[,c("crAssphage","Gubaphage")])
chisq.test(chi.test.data[,c("High-quality","Gubaphage")])
mt2<-genome.info %>% 
  mutate(MTase=factor(MTase,levels = c("Y","N"))) %>% 
  mutate(Type=factor(Type,levels = c("High-quality","crAssphage","Gubaphage"))) %>% 
  
  ggplot()+
  geom_bar(aes(x=Type,y=..count..,fill=MTase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+
labs(y="Percentage",title = "MTase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))

chi.test.data<-aggregate(genome.info$REase,by=list(genome.info$Type,genome.info$REase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chisq.test(chi.test.data[,c("High-quality","crAssphage")])
chisq.test(chi.test.data[,c("crAssphage","Gubaphage")])
chisq.test(chi.test.data[,c("High-quality","Gubaphage")])
re2<-genome.info %>% 
  mutate(REase=factor(REase,levels = c("Y","N"))) %>% 
  mutate(Type=factor(Type,levels = c("High-quality","crAssphage","Gubaphage"))) %>% 
  ggplot()+
  geom_bar(aes(x=Type,y=..count..,fill=REase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+
  labs(y="",title = "REase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))
(mt+re)/(mt2+re2)
#Figure2D----
chi.test.data<-aggregate(GPD.meta$MTase,by=list(GPD.meta$Lifestyle,GPD.meta$MTase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chisq.test(chi.test.data[,c(2,3)])
chisq.test(chi.test.data[,c(3,4)])
chisq.test(chi.test.data[,c(4,5)])
mt<-GPD.meta %>% 
  subset(Anno=="High-quality") %>%
  mutate(MTase=factor(MTase,levels = c("Y","N"))) %>% 
  ggplot()+
  geom_bar(aes(x=Lifestyle,y=..count..,fill=MTase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+

labs(y="Percentage",title = "GPD MTase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))

chi.test.data<-aggregate(GPD.meta$REase,by=list(GPD.meta$Lifestyle,GPD.meta$REase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chisq.test(chi.test.data[,c(2,5)])
chisq.test(chi.test.data[,c(2,3)])
chisq.test(chi.test.data[,c(3,4)])
chisq.test(chi.test.data[,c(4,5)])
re<-GPD.meta %>% 
  subset(Anno=="High-quality") %>%
  
  mutate(REase=factor(REase,levels = c("Y","N"))) %>% 
  ggplot()+
  geom_bar(aes(x=Lifestyle,y=..count..,fill=REase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+

labs(y="",title = "GPD REase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))
# mt+re
chi.test.data<-aggregate(GPD.meta$MTase,by=list(GPD.meta$Type,GPD.meta$MTase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chi.test.data[is.na(chi.test.data)]<-0
chisq.test(chi.test.data[,c("High-quality","crAssphage")])
chisq.test(chi.test.data[,c("crAssphage","Gubaphage")])
chisq.test(chi.test.data[,c("High-quality","Gubaphage")])
mt2<-GPD.meta %>% 
  mutate(MTase=factor(MTase,levels = c("Y","N"))) %>% 
  mutate(Type=factor(Type,levels = c("High-quality","crAssphage","Gubaphage"))) %>% 
  subset(Anno=="High-quality") %>%
  
  ggplot()+
  geom_bar(aes(x=Type,y=..count..,fill=MTase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+

labs(y="Percentage",title = "GPD MTase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))

chi.test.data<-aggregate(GPD.meta$REase,by=list(GPD.meta$Type,GPD.meta$REase),FUN=length)
chi.test.data<-dcast(chi.test.data,Group.2~Group.1)
chisq.test(chi.test.data[,c("High-quality","crAssphage")])
chisq.test(chi.test.data[,c("crAssphage","Gubaphage")])
chisq.test(chi.test.data[,c("High-quality","Gubaphage")])
re2<-GPD.meta %>% 
  subset(Anno=="High-quality") %>%
  
  mutate(REase=factor(REase,levels = c("Y","N"))) %>% 
  mutate(Type=factor(Type,levels = c("High-quality","crAssphage","Gubaphage"))) %>% 
  ggplot()+
  geom_bar(aes(x=Type,y=..count..,fill=REase),position = "fill",width=0.7)+
  scale_fill_manual(values = c("#3B4992","#89A2E5","#A76CA8","#DD5E5E"))+
  theme_bw()+
labs(y="",title = "GPD REase")+
  guides(fill="none")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.10),breaks=seq(0,1,0.25))+
  theme(axis.text.x = element_text(angle = 15,hjust = 1))

(mt+re)/(mt2+re2)
#Figure3B----
MTase.prophage.blast<-MTase.prophage.blast%>% group_by(pt1) %>% 
  filter(similarity==max(similarity)) %>% ungroup()
total<-1281
o90<-subset(MTase.prophage.blast,similarity>=90)$pt1 %>% unique%>% length
o70<-subset(MTase.prophage.blast,similarity>=70 & similarity<90)$pt1 %>% unique%>% length
o50<-subset(MTase.prophage.blast,similarity>=50 & similarity<70)$pt1 %>% unique%>% length
o30<-subset(MTase.prophage.blast,similarity>=30 & similarity<50)$pt1 %>% unique%>% length
simi2bac<-data.frame(over90=o90/total,
                     over70=o70/total,
                     over50=o50/total,
                     over30=o30/total,
                     
                     group="MTase")

RE.prophage.blast<-RE.prophage.blast%>% group_by(pt1) %>% 
  filter(similarity==max(similarity)) %>% ungroup()
total<-623
o90<-subset(RE.prophage.blast,similarity>=90)$pt1 %>% unique%>% length
o70<-subset(RE.prophage.blast,similarity>=70 & similarity<90)$pt1 %>% unique%>% length
o50<-subset(RE.prophage.blast,similarity>=50 & similarity<70)$pt1 %>% unique%>% length
o30<-subset(RE.prophage.blast,similarity>=30 & similarity<50)$pt1 %>% unique%>% length
simi2bac<-rbind(simi2bac,data.frame(over90=o90/total,
                                    over70=o70/total,
                                    over50=o50/total,
                                    over30=o30/total,
                                    group="REase"))
simi2bac$less30<-1-rowSums(simi2bac[,1:4])
simi2bac<-melt(simi2bac)
simi2bac$variable<-gsub("over90",">=90",simi2bac$variable)
simi2bac$variable<-gsub("over70","70~90",simi2bac$variable)
simi2bac$variable<-gsub("over50","50~70",simi2bac$variable)
simi2bac$variable<-gsub("over30","30~50",simi2bac$variable)
simi2bac$variable<-gsub("less30","<30",simi2bac$variable)

simi2bac$variable<-factor(simi2bac$variable,levels=rev(c(">=90","70~90",
                                                         "50~70","30~50","<30")))
simi2bac = simi2bac[order(simi2bac$variable, decreasing = FALSE),]
myLabel = as.vector(subset(simi2bac,group=="MTase")$variable)   
myLabel = paste(myLabel, "(", round(subset(simi2bac,group=="MTase")$value * 100, 2), "%)", sep = "")   
p1<-subset(simi2bac,group=="MTase") %>% 
  mutate(variable=fct_reorder(variable,value)) %>% 
  ggplot( aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  
  theme_bw()+
  labs(y="",fill="Similarities to B-MTase")+
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) + 
  facet_grid(~"MVP V-MTase")+ 
  scale_fill_brewer(labels = myLabel) +
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))


simi2bac$variable<-factor(simi2bac$variable,levels=rev(c(">=90","70~90",
                                                         "50~70","30~50","<30")))
simi2bac = simi2bac[order(simi2bac$variable, decreasing = FALSE),]
myLabel2 = as.vector(subset(simi2bac,group=="REase")$variable)   
myLabel2 = paste(myLabel2, "(", round(subset(simi2bac,group=="REase")$value * 100, 2), "%)", sep = "")   
p2<-subset(simi2bac,group=="REase") %>% 
  mutate(variable=factor(variable,levels=rev(c(">=90","70~90",
                                               "50~70","30~50","<30")))) %>% 
  ggplot( aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  
  theme_bw()+
  labs(y="",fill="Similarities to B-REase")+
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) + 
  facet_grid(~"MVP V-REase")+ 
  scale_fill_brewer(palette = "YlOrRd" ,labels = myLabel2) +
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))
p1/p2
#Figure3C----
load("Data/mvp.blast.RData")
#species--
RE.prophage.blast$"WithHost"<-"N"
RE.prophage.blast$"WithHost"[which(RE.prophage.blast$s1==RE.prophage.blast$s2 )]<-"Y"
MTase.prophage.blast$"WithHost"<-"N"
MTase.prophage.blast$"WithHost"[which(MTase.prophage.blast$s1==MTase.prophage.blast$s2 )]<-"Y"
RE.prophage.blast$"Group"<-"REase"
MTase.prophage.blast$"Group"<-"MTase"
prophage.blast<-rbind(MTase.prophage.blast,RE.prophage.blast) %>% unique
multi<-(paste(prophage.blast$ge1,prophage.blast$s2,sep="__") %>% table %>% data.frame %>% subset(Freq>1))[1]
colnames(multi)<-"l"
prophage.blast.multi<- subset(prophage.blast,paste(ge1,s1,sep="__") %in% as.character(multi$l))
# prophage.blast.multi<-prophage.blast.multi %>% group_by(pt1) %>%
#   filter(similarity==max(similarity)) %>% ungroup()
plot.roc(subset(prophage.blast.multi,Group=="REase" & Quality=="High" )$WithHost, 
         subset(prophage.blast.multi,Group=="REase" & Quality=="High" )$similarity,
         percent=TRUE, thresholds="best", # select the (best) threshold
         print.thres="best",
         ci=TRUE,
         print.auc=TRUE, #display pAUC value on the plot with following options:
         auc.polygon=TRUE, auc.polygon.col="#BF5454", # show pAUC as a polygon
         max.auc.polygon=TRUE, max.auc.polygon.col="#BF545422", # also show the 100% polygon
         main="REase ROC plot(species level)")
plot.roc(subset(prophage.blast.multi,Group=="MTase" & Quality=="High" )$WithHost, 
         subset(prophage.blast.multi,Group=="MTase" & Quality=="High" )$similarity,
         percent=TRUE, thresholds="best", # select the (best) threshold
         print.thres="best",
         ci=TRUE,
         print.auc=TRUE, #display pAUC value on the plot with following options:
         auc.polygon=TRUE, auc.polygon.col="#1c61b6", # show pAUC as a polygon
         max.auc.polygon=TRUE, max.auc.polygon.col="#1c61b622", # also show the 100% polygon
         main="MTase ROC plot(species level)")
#Genus--
RE.prophage.blast$"WithHost"<-"N"
RE.prophage.blast$"WithHost"[which(RE.prophage.blast$s1==RE.prophage.blast$s2 |RE.prophage.blast$g1==RE.prophage.blast$g2)]<-"Y"
MTase.prophage.blast$"WithHost"<-"N"
MTase.prophage.blast$"WithHost"[which(MTase.prophage.blast$s1==MTase.prophage.blast$s2 |MTase.prophage.blast$g1==MTase.prophage.blast$g2)]<-"Y"
RE.prophage.blast$"Group"<-"REase"
MTase.prophage.blast$"Group"<-"MTase"
prophage.blast<-rbind(MTase.prophage.blast,RE.prophage.blast) %>% unique
multi<-(paste(prophage.blast$ge1,prophage.blast$s2,sep="__") %>% table %>% data.frame %>% subset(Freq>1))[1]
colnames(multi)<-"l"
prophage.blast.multi<- subset(prophage.blast,paste(ge1,s1,sep="__") %in% as.character(multi$l))
prophage.blast.multi<-prophage.blast.multi %>% group_by(pt1) %>% 
  filter(similarity==max(similarity)) %>% ungroup()
plot.roc(subset(prophage.blast.multi,Group=="REase" & Quality=="High" )$WithHost, 
         subset(prophage.blast.multi,Group=="REase" & Quality=="High" )$similarity,
         percent=TRUE, thresholds="best", # select the (best) threshold
         print.thres="best",
         ci=TRUE,
         print.auc=TRUE, #display pAUC value on the plot with following options:
         auc.polygon=TRUE, auc.polygon.col="#BF5454", # show pAUC as a polygon
         max.auc.polygon=TRUE, max.auc.polygon.col="#BF545422", # also show the 100% polygon
         main="REase ROC plot(genus level)")
plot.roc(subset(prophage.blast.multi,Group=="MTase" & Quality=="High" )$WithHost, 
         subset(prophage.blast.multi,Group=="MTase" & Quality=="High" )$similarity,
         percent=TRUE, thresholds="best", # select the (best) threshold
         print.thres="best",
         ci=TRUE,
         print.auc=TRUE, #display pAUC value on the plot with following options:
         auc.polygon=TRUE, auc.polygon.col="#1c61b6", # show pAUC as a polygon
         max.auc.polygon=TRUE, max.auc.polygon.col="#1c61b622", # also show the 100% polygon
         main="MTase ROC plot(genus level)")

#Figure3D----
RE.prophage.blast$"WithHost"<-"N"
RE.prophage.blast$"WithHost"[which(RE.prophage.blast$s1==RE.prophage.blast$s2 )]<-"Y"
MTase.prophage.blast$"WithHost"<-"N"
MTase.prophage.blast$"WithHost"[which(MTase.prophage.blast$s1==MTase.prophage.blast$s2 )]<-"Y"
RE.prophage.blast$"Group"<-"REase"
MTase.prophage.blast$"Group"<-"MTase"
prophage.blast<-rbind(MTase.prophage.blast,RE.prophage.blast) %>% subset(similarity>=90) %>% unique
multi<-(paste(prophage.blast$ge1,prophage.blast$s2,sep="__") %>% table %>% data.frame %>% subset(Freq>1))[1]
colnames(multi)<-"l"
prophage.blast.multi<- subset(prophage.blast,paste(ge1,s1,sep="__") %in% as.character(multi$l))
prophage.blast.multi<-prophage.blast.multi %>% group_by(pt1) %>% 
  filter(similarity==max(similarity)) %>% ungroup()

p1<-data.frame(WithHost=c("N","Y"),
               count=c((subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="N")$ge1 %>% unique %>% length),
                       (subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="Y")$ge1 %>% unique %>% length)),
               GeneType=c("MTase/REase","MTase/REase")) %>% 
  ggplot(aes(x=GeneType,y=count,fill=WithHost))+
  geom_bar(position="fill", stat = "identity",alpha=0.9,width = 0.5) +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values=c("#ebcdc5","#910c29"))+
  theme_bw()+
  facet_grid(~"Species level")+
  labs(y="Specificity",fill="")+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))
data.frame(WithHost=c("N","Y"),
           count=c((subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="N")$ge1 %>% unique %>% length),
                   (subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="Y")$ge1 %>% unique %>% length)),
           GeneType=c("MTase/REase","MTase/REase"))

RE.prophage.blast$"WithHost"<-"N"
RE.prophage.blast$"WithHost"[which(RE.prophage.blast$s1==RE.prophage.blast$s2 |RE.prophage.blast$g1==RE.prophage.blast$g2)]<-"Y"
MTase.prophage.blast$"WithHost"<-"N"
MTase.prophage.blast$"WithHost"[which(MTase.prophage.blast$s1==MTase.prophage.blast$s2 |MTase.prophage.blast$g1==MTase.prophage.blast$g2)]<-"Y"
RE.prophage.blast$"Group"<-"REase"
MTase.prophage.blast$"Group"<-"MTase"
prophage.blast<-rbind(MTase.prophage.blast,RE.prophage.blast) %>% subset(similarity>=90) %>% unique
multi<-(paste(prophage.blast$ge1,prophage.blast$s2,sep="__") %>% table %>% data.frame %>% subset(Freq>1))[1]
colnames(multi)<-"l"
prophage.blast.multi<- subset(prophage.blast,paste(ge1,s1,sep="__") %in% as.character(multi$l))
prophage.blast.multi<-prophage.blast.multi %>% group_by(pt1) %>% 
  filter(similarity==max(similarity)) %>% ungroup()

p2<-data.frame(WithHost=c("N","Y"),
               count=c((subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="N")$ge1 %>% unique %>% length),
                       (subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="Y")$ge1 %>% unique %>% length)),
               GeneType=c("MTase/REase","MTase/REase")) %>% 
  ggplot(aes(x=GeneType,y=count,fill=WithHost))+
  geom_bar(position="fill", stat = "identity",alpha=0.9,width = 0.5) +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values=c("#ebcdc5","#910c29"))+
  theme_bw()+
  facet_grid(~"Genus level")+
  labs(y="Specificity",fill="")+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))
data.frame(WithHost=c("N","Y"),
           count=c((subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="N")$ge1 %>% unique %>% length),
                   (subset(prophage.blast.multi,Quality=="High" & similarity>=90 &WithHost=="Y")$ge1 %>% unique %>% length)),
           GeneType=c("MTase/REase","MTase/REase"))
p1+p2+plot_layout(guides = "collect")

#Figure3E----
ases.blast.filtered<-ases.blast%>% subset(similarity>=90) %>% unique 
multi<-(paste(ases.blast.filtered$g1,ases.blast.filtered$Species,sep="__") %>% table %>% data.frame %>% subset(Freq>1))[1]
colnames(multi)<-"l"

ases.blast.filtered.multi<- subset(ases.blast.filtered,paste(g1,Species,sep="__") %in% as.character(multi$l))
ases.blast.filtered<-ases.blast.filtered.multi
ases.blast.filtered$"Type"<-"High-quality"
ases.blast.filtered$"Type"[which(ases.blast.filtered$g1 %in%crass.list$V1 )]<-"crAssphage"
ases.blast.filtered$"Type"[which(ases.blast.filtered$g1 %in%guba.list$V1 )]<-"Gubaphage"
num<-(ases.blast.filtered%>% select(g1) %>% unique)$g1 %>% length
num1<-(ases.blast.filtered %>% subset(Group=="MTase" )%>% select(g1) %>% unique)$g1 %>% length
num2<-(ases.blast.filtered %>% subset(Group=="REase" )%>% select(g1) %>% unique)$g1 %>% length
num3<-subset(ases.blast.filtered,Group=="MTase" & g1 %in% subset(ases.blast.filtered,Group=="REase")$g1)$g1 %>% unique %>% length

data=data.frame(variable=c("MTase","Both","REase","Not available"),
                value=c(num1-num3,num3,num2-num3,9665-num))
myLabel = as.vector(data$variable)   
myLabel = paste(myLabel, " (", round(data$value/sum(data$value) * 100, 2), "%)", sep = "")   
data %>% 
  mutate(variable=factor(variable,levels = c("MTase","Both","REase","Not available"))) %>% 
  ggplot( aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.1,color="white",size=1.5) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme_void()+
  labs(y="",fill="")+
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position = "bottom") + 
  guides(fill = guide_legend( ncol = 1),color="none")+
  scale_fill_manual(values=c("#42777a","#432c39","#d74c35","#efd9b5"),labels = myLabel) +
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))+
  theme(legend.text=element_text(size=14))

#Figure3F----
#Host range MTases--
color<-c("#648E8C","#A6BEB7","#EBE25B","#FF95A8FF","#FF6348FF","#CF5845","#A71F28") %>% rev()
host2v<-subset(ases.blast.filtered,Group=="MTase")[,c("g1","Kingdom","Phylum","Class","Order","Family","Genus","Species")] 
colnames(host2v)[1]<-"Virus"
host2v.LCA<-calculate_LCA(host2v)
HostRange<-table(host2v.LCA$level) %>% data.frame()
HostRange$"Level"<-"MTase"
names(HostRange)<-c("Range","Count","Level")
myLabel = as.vector(HostRange$Range)   
myLabel = paste(myLabel, "(", round(HostRange$Count/ sum(HostRange$Count) * 100, 2), "%)", sep = "")   
HostRange$"label"<-myLabel
HostRange<-HostRange %>%
  mutate(Range=factor(Range,levels = c("Kingdom","Phylum","Class","Order",
                                       "Family","Genus","Species"))) 
mylable<-arrange(HostRange,Range)$label
HostRange<-HostRange %>%
  mutate(label=factor(label,levels=mylable))
p1<-HostRange %>%
  ggplot(aes(x=1,y = Count, fill = label)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "",fill="Host Range") + 
  theme(axis.ticks = element_blank()) + 
  theme_bw()+
  scale_fill_brewer(palette = "RdYlBu")+

  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  facet_grid(~Level)+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))
#Host range REases--

host2v<-subset(ases.blast.filtered,Group=="REase")[,c("g1","Kingdom","Phylum","Class","Order","Family","Genus","Species")] 
colnames(host2v)[1]<-"Virus"
host2v.LCA<-calculate_LCA(host2v)
HostRange<-table(host2v.LCA$level) %>% data.frame()
HostRange$"Level"<-"REase"
names(HostRange)<-c("Range","Count","Level")
myLabel = as.vector(HostRange$Range)   
myLabel = paste(myLabel, "(", round(HostRange$Count/ sum(HostRange$Count) * 100, 2), "%)", sep = "")   
HostRange$"label"<-myLabel
HostRange<-HostRange %>%
  mutate(Range=factor(Range,levels = c("Kingdom","Phylum","Class","Order",
                                       "Family","Genus","Species"))) 
mylable<-arrange(HostRange,Range)$label
HostRange<-HostRange %>%
  mutate(label=factor(label,levels=mylable))
p2<-HostRange %>%
  ggplot(aes(x=1,y = Count, fill = label)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "",fill="Host Range") + 
  theme(axis.ticks = element_blank()) + 
  theme_bw()+
  scale_fill_brewer(palette = "RdYlBu")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))+
  facet_grid(~Level)
#Host range All--
host2v<-ases.blast.filtered[,c("g1","Kingdom","Phylum","Class","Order","Family","Genus","Species")] 
colnames(host2v)[1]<-"Virus"
host2v.LCA<-calculate_LCA(host2v)
HostRange<-table(host2v.LCA$level) %>% data.frame()
HostRange$"Level"<-"Combined"
names(HostRange)<-c("Range","Count","Level")
myLabel = as.vector(HostRange$Range)   
myLabel = paste(myLabel, "(", round(HostRange$Count/ sum(HostRange$Count) * 100, 2), "%)", sep = "")   
HostRange$"label"<-myLabel
HostRange<-HostRange %>%
  mutate(Range=factor(Range,levels = c("Kingdom","Phylum","Class","Order",
                                       "Family","Genus","Species"))) 
mylable<-arrange(HostRange,Range)$label
HostRange<-HostRange %>%
  mutate(label=factor(label,levels=mylable))
p3<-HostRange %>%
  ggplot(aes(x=1,y = Count, fill = label)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "",fill="Host Range") + 
  theme(axis.ticks = element_blank()) + 
  theme_bw()+
  scale_fill_brewer(palette = "RdYlBu")+

  theme(strip.background.x = element_rect(fill = "white"),strip.background.y = element_rect(fill = "white"))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  facet_grid(~Level)
#Host range plot--
p1/p2/p3
