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
cluster_color = list(
   veg = "darkorange1",
   early = "deepskyblue",
   inter = "chartreuse3",
   late = "red",
   very_late = "deeppink2"
 )

Culster_color <- function(data_tab, list = rnai_list, cluster_list = cluster, color_list = cluster_color){
  color = c()
  
  rnai = unique(str_remove_all(str_split_fixed(colnames(data_tab), "_T", n=2)[,1], "_Veg"))

  for (r in rnai){
    for (clust in cluster_list[[r]]){
      if (clust == "VEG"){
        color = c(color, color_list[["veg"]])
      } else if(clust == "EARLY" ){
        color = c(color, color_list[["early"]])
      } else if(clust == "INTER" ){
        color = c(color, color_list[["inter"]])
      } else if(clust == "LATE" ){
        color = c(color, color_list[["late"]])
      }
    }
  }
  return(color)
}

Culster_color_info <- function(data_tab, infodata , color_list = cluster_color){
  
  if(is.element(FALSE, colnames(data_tab)==rownames(infodata))){
    infodata = infodata[colnames(data_tab),]
  }
  cluster = infodata$Cluster
  color = c()
    for (clust in cluster){
      if (clust == "VEG"){
        color = c(color, color_list[["veg"]])
      } else if(clust == "EARLY" ){
        color = c(color, color_list[["early"]])
      } else if(clust == "INTER" ){
        color = c(color, color_list[["inter"]])
      } else if(clust == "LATE" ){
        color = c(color, color_list[["late"]])
      }
    }
  names(color) = rownames(infodata)

  return(color)
}

method_color = list(
  HiSeq = "chartreuse4",
  NextSeq = "blue4",
  both = "mediumturquoise")


Batch_color <- function(data_tab, list= rnai_list, cluster_list = cluster, color_list = method_color){
  rnai = unique(str_remove_all(str_split_fixed(colnames(data_tab), "_T", n=2)[,1], "_Veg"))
  
  seq_color = c()
  
  if(any(grepl("bis",rnai))){
    for (j in rnai){
      if (is.element(j,c("ICL7bis", "EZL1bis", "XRCC4", "CTIP", "ND7_C", "ND7_X"))){
        seq_color=c(seq_color,rep(color_list[["NextSeq"]],length(cluster_list[[j]])))
      }else if (is.element(j, c("ICL7", "EZL1", "ND7_K", "KU80c", "PGM"))){
        seq_color=c(seq_color,rep(color_list[["HiSeq"]], length(cluster_list[[j]])))
      }
    }
  }else{
    for (j in rnai){
      if (is.element(j,c("ICL7", "EZL1"))){
        if(length(cluster_list[["ICL7"]])==length(cluster_list[["ICL7bis"]])){
          seq_color=c(seq_color,rep(color_list[["NextSeq"]],length(cluster_list$ICL7bis)))
        }else{
          seq_color=c(seq_color,
                      method_color$both, method_color$both, method_color$HiSeq,
                      method_color$both, method_color$HiSeq, method_color$both, method_color$HiSeq)
        }
        
      }else if (is.element(j,c("XRCC4", "CTIP", "ND7_C", "ND7_X"))){
        seq_color=c(seq_color,rep(color_list[["NextSeq"]],length(cluster_list[[j]])))
      }else if (is.element(j, c("ND7_K", "KU80c", "PGM"))){
        seq_color=c(seq_color,rep(color_list[["HiSeq"]], length(cluster_list[[j]])))
      }
    }
  }
  
  return(seq_color)
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
rm(selection,i,name)


