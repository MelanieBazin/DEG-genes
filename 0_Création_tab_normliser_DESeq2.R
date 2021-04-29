source("0_Visualisatin  des données_fonction.R", encoding = "UTF-8")
source("2_Fonction.R")

# RNAi à analyser ensemble
tout = sub(paste0("_expression_table_RPKM.tab"),"",list.files("./DATA/RPKM/"))
rnai_list = list(
  tout = tout,
  sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
  sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
  controles = tout[which(is.element(tout,c("ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
  XRCC4ctrl2020 = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","CTIP_CTRL")))],
  XRCC4tousctrl = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","ND7", "ICL7", "CTIP_CTRL")))],
  XRCC4xseq2014 = tout[which(is.element(tout,c("XRCC4","XRCC4_CTRL","ICL7","KU80c","ND7","PGM" )))],
  CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
  CTIPctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL","XRCC4_CTRL")))],
  CTIPtousctrl = tout[which(is.element(tout,c("CTIP","ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))]
)
i = names(rnai_list)[3]

for (i in names(rnai_list)){
  #Ouverture des fichiers et création de l'objet countdata
  countdata = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
  row.names(countdata)=countdata$ID
  countdata=countdata[,-1]
  countdata =  as.matrix(countdata)
  
  # Création du tableau avec les info des colonnes
  infodata=CreatInfoData2(conditions = rnai_list[[i]])
  
  # Mise en forme des données
  library(DESeq2)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ RNAi + Timing)
  
  
  # Analyse DESeq2
  deseq = DESeq(deseq)
  
  # Graphique du paramètre de dispersion
  png(paste0("./Graph/DESeq2/",i,"_dipression_DESeq2.png"))
    plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()
  
  # Récupération des données de comptage normalisées
  tab=counts(deseq,normalized=T)
  write.table(tab,paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_normalisation_DESeq2.tab"), row.names = T, sep="\t")
  write.table(infodata,paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_infodata_DESeq2.tab"), row.names = T, sep="\t")
  
}

# normacount = read.table(paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_normalisation_DESeq2.tab"), header = T, sep="\t")
# info = read.table(paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_infodata_DESeq2.tab"), header = T, sep="\t")

