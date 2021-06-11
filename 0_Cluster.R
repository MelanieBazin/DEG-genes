options(stringsAsFactors = FALSE)
library(stringr)
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
rownames(annotation)=annotation$ID

#### Définition des ARNi à analyser ensembles ####
tout = sub(".tab","",list.files("./DATA/EXPRESSION/"))

rnai_list = list(
  sequencage_2014 = tout[which(is.element(tout,c("KU80c","ND7_K","PGM", "ICL7" )))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","ND7_X")))],
  CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","ND7_C", "ND7_X")))]
)
rm(tout)

#### Définition des cluster à grouper ensemble ####
cluster = list(
  ICL7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  ND7_K = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  PGM = c(rep("VEG",1),rep("INTER",3),rep("LATE",3)),
  KU80c = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  
  ND7_C = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",2)),
  CTIP = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",1)),
  ND7_X = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
  XRCC4 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2))
)


#### Définition des couleur à attribuer pour les différents RNAi ###"
veg_color = "darkorange1"
early_color = "deepskyblue"
inter_color = "chartreuse3"
late_color = "red"
very_late_color = "deeppink2"

cluster_color = list()
for(j in names(cluster)){
  vec = cluster[[j]]
  color = c()
  for (i in 1:length(vec)){
    if (vec[i] == "VEG"){
      color=c(color,veg_color)
    }else if (vec[i] == "EARLY"){
      color=c(color,early_color)
    }else if (vec[i] == "INTER"){
      color=c(color,inter_color)
    }else if (vec[i] == "LATE"){
      color=c(color,late_color)
    }else if (vec[i] == "VERY_LATE"){
      color=c(color,very_late_color)
    }
  }
  cluster_color[[j]] = color
}
rm(veg_color, early_color, inter_color, late_color, very_late_color,i,j)


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


