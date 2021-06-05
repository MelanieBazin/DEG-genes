options(stringsAsFactors = FALSE)
library(stringr)

#### Définition des ARNi à analyser ensembles ####
tout = sub(".tab","",list.files("./DATA/EXPRESSION/"))

rnai_list = list(
  tout = tout,
  sequencage_2014 = tout[which(is.element(tout,c("KU80c","ND7_K","PGM", "ICL7" )))],
  sequencage_2014bis = tout[which(is.element(tout,c("KU80c","ND7_K","PGM" )))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","ND7_X")))],
  CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","ND7_C", "ND7_X")))]
)
rm(tout)

#### Définition des cluster à grouper ensemble ####
cluster = list(
  ICL7 = c(rep("EARLY",1),rep("INTER",2),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  ND7_K = c(rep("EARLY",1),rep("INTER",2),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  PGM = c(rep("INTER",3),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  KU80c = c(rep("EARLY",1),rep("INTER",2),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  
  ND7_C = c(rep("EARLY",1),rep("INTER",3),rep("LATE",1),rep("VEG",1)),
  CTIP = c(rep("EARLY",1),rep("INTER",2),rep("LATE",1),rep("VEG",1)),
  ND7_X = c(rep("EARLY",1),rep("INTER",2),rep("LATE",1),rep("VEG",1)),
  XRCC4 = c(rep("EARLY",1),rep("INTER",2),rep("LATE",1),rep("VEG",1))
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
rm(veg_color, early_color, inter_color, late_color, very_late_color)


#### Definition de l'ordre des colonnes #####
tabs = list.files("./DATA/EXPRESSION")
timing_list = as.list(tabs)
for (k in 1:length(tabs)){
  table = read.table(paste0("./DATA/EXPRESSION/",tabs[k]), header = T, row.names = 1)
  names(timing_list)[k] = gsub(".tab", "", tabs[k])
  timing_list[[k]] = colnames(table)
}
rm(table, tabs)
