source("0_Cluster.R")
source("3_Visualisation_des_donnees_fonction.R")
library(sva)
library(DESeq2)
tout = sub(".tab","",list.files("./DATA/EXPRESSION/"))
rnai_list = list(
  tout = tout,
  sequencage_2014 = tout[which(is.element(tout,c("KU80c","ND7_K","PGM", "ICL7" )))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","ND7_X")))],
  CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","ND7_C", "ND7_X")))]
)

seq_2014 = "chartreuse4"
seq_2020 = "blue4"

seq_color = list()
for(j in names(cluster)){
  vec = cluster[[j]]
  color = c()
    if (is.element(j, c("KU80c","ND7_K","PGM", "ICL7" ))){
      color=c(color,rep(seq_2014,length(vec)))
    }else if (is.element(j, c("XRCC4","ND7_X","CTIP","ND7_C"))){
      color=c(color,rep(seq_2020, length(vec)))
    }
  
  seq_color[[j]] = color
}

condition = "tout"
#### ACP sans correction Batch ####
countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[condition]])

infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)
countdata = as.matrix(countdata)

deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition
                               # design   = ~ Feeding + Cluster
)

deseq = DESeq(deseq)
data_tab = counts(deseq,normalized=T)
# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(j, color)]=seq_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = "./Analyse/ACP_tout/",
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

#### Refaire ACP avec la correction Batch ####
countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[condition]])
infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)
countdata = as.matrix(countdata)
countdata = ComBat_seq(countdata, batch = infodata$Batch, group = infodata$Cluster)



deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition
                               # design   = ~ Feeding + Cluster
)

deseq = DESeq(deseq)
data_tab = counts(deseq,normalized=T)

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = "./Analyse/ACP_tout/Batch_color",
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(j, color)]=cluster_color[[j]]
}
