options(stringsAsFactors = FALSE)
library(stringr)
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
rownames(annotation)=annotation$ID

#### Définition des ARNi à analyser ensembles ####
tout = sub(".tab","",list.files("./DATA/EXPRESSION/"))

rnai_list = list(
  HiSeqvsNextSeq = tout[which(is.element(tout,c("EZL1bis","ICL7bis","EZL1", "ICL7" )))],
  analyseDE = tout[which(!is.element(tout,c("EZL1bis","EZL1" )))]
)
rm(tout)

#### Définition des cluster à grouper ensemble ####
cluster = list(
  ICL7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  ICL7bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
  EZL1 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  EZL1bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
  
  ND7_K = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  PGM = c(rep("VEG",1),rep("INTER",3),rep("LATE",3)),
  KU80c = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  
  ND7_C = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",2)),
  CTIP = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",1)),
  ND7_X = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
  XRCC4 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2))
)


#### Définition des couleur à attribuer pour les différents RNAi ###"
 clust_color = list(
   veg = "darkorange1",
   early = "deepskyblue",
   inter = "chartreuse3",
   late = "red",
   very_late = "deeppink2"
 )

Culster_color <- function(data_tab, infodata, clust_color){
  clus = c()
  for (c in colnames(data_tab)){
    pos = which(infodata$Names == c)
    clus = c(clus, infodata$Cluster[pos])
  }
  
  for (c in names(clust_color)){
    pos = grep(c, clus, ignore.case = T)
    clus[pos] = clust_color[[c]]
  }
  return(clus)
}

batch_color = list(
  seq_2014 = "chartreuse4",
  seq_2020 = "blue4",
  both = "mediumturquoise")


Batch_color <- function(data_tab, infodata, batch_color){
  clus = c()
  for (c in colnames(data_tab)){
    pos = which(infodata$Names == c)
    if(grepl(",",infodata$runsCollapsed[pos])){
      clus = c(clus, "both")
    }else{
      clus = c(clus, infodata$Batch[pos])
    }
  }
  
  for (c in names(batch_color)){
    pos = grep(c, clus, ignore.case = T)
    clus[pos] = batch_color[[c]]
  }
  return(clus)
}

#### Definition de l'ordre des colonnes #####
tabs = list.files("./DATA/EXPRESSION")
timing_list = as.list(tabs)
for (k in 1:length(tabs)){
  table = read.table(paste0("./DATA/EXPRESSION/",tabs[k]), header = T, row.names = 1)
  names(timing_list)[k] = gsub(".tab", "", tabs[k])
  timing_list[[k]] = colnames(table)
}
rm(table, tabs,k)

####### Sélection de gènes #####
selection = c("Ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")
selection = sort(selection)

select_ID =c()
name = c()
for( i in 1:length(selection)){
  select_ID = c(select_ID,annotation$ID[grep(selection[i],annotation$NAME, ignore.case = T)]) 
  name = c(name,annotation$NAME[grep(selection[i],annotation$NAME, ignore.case = T)]) 
}
names(select_ID)=name
# select_annotation = annotation[select_ID,]
rm(selection,i,name)


