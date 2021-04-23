options(stringsAsFactors = FALSE)

annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

type = "EXPRESSION" #3 possibiltées : "RPM","RPKM", "EXPRESSION"

condition = "CTIP"

#### Ouverture des fichiers ####
files = list.files(path = paste0("./DATA/",type))
files = files[grep(condition,files)]

l_tab = vector("list",length(files))
for (i in 1:length(files)){
  
  tab = read.table(paste0("./DATA/", type,"/",files[i]),h = T,sep = "\t")
  tab = merge(annotation,tab, x.by = "ID", y.by = "ID")
  
  l_tab[[i]] = tab
  
}

control = l_tab[[grep("CTRL",files)]]
rnai = l_tab[[setdiff(1:length(files),grep("CTRL",files))]]

### Mise en forme des tableau et données pour graphiques ########
if (condition == "CTIP") {
  timing_ctrl = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",2))
  timing_rnai = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",1))
  
  control = control[,c(1:6,12,7:11)]
  rnai = rnai[,c(1:6,11,7:10)]
  FC = -0.1
  
}
if (condition == "XRCC4"){
  timing_ctrl = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2))
  timing_rnai = timing_ctrl
  
  control = control[,c(1:6,11,7:10)]
  rnai = rnai[,c(1:6,11,7:10)]
  FC = NULL
}

####### Sélection de gènes #####
selection = c("Ku","ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")
selection = sort(selection)

select_ctrl =c()
select_rnai = c()

for( i in 1:length(selection)){
  select_ctrl = c(select_ctrl,grep(selection[i],control$SYNONYMS)) 
  select_rnai = c(select_rnai,grep(selection[i],rnai$SYNONYMS))
  }
select_ctrl = control[select_ctrl,]
select_rnai = rnai[select_rnai,]

####### Mise en forme des données pour DESeq2 ##############
select_ctrl = select_ctrl[,c(1,7:ncol(select_ctrl))]
select_rnai = select_rnai[,c(1,7:ncol(select_rnai))]
colnames(select_ctrl) = paste(colnames(select_ctrl),"CTRL", sep = "_" )
colnames(select_ctrl)[1] = "ID"
colnames(select_rnai) = paste(colnames(select_rnai),"RNAi", sep = "_" )
colnames(select_rnai)[1] = "ID"

countdata = merge(select_ctrl,select_rnai, x.by = "ID", y.by = "ID")
row.names(countdata)=countdata$ID[]
countdata = as.matrix(countdata[,-1])

infodata = matrix(NA,c(ncol(countdata),01))
infodata[,1] = c(rep("controle",ncol(select_ctrl)-1),rep(paste("RNAi",condition),ncol(select_rnai)-1))
infodata = cbind(infodata, c(timing_ctrl,timing_rnai))

infodata = cbind(colnames(countdata), infodata)
colnames(infodata) = c("noms", "condition", "timing")
row.names(infodata)=infodata[,1]

##### Chargement des librairies pour l'analyse #############
library(FactoMineR)
library(factoextra)
library(DESeq2)

##### Premier représentataion des données ###########
  
# Filtration des gènes avec trop peu de compatage (seuil arbitraire)
countdata = countdata[rowSums(countdata) > 50,]


## Nomralisation via DESeq2
# Mise en forme des données
coldata = data.frame(infodata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                                      colData  = coldata,
                                      design   = ~ condition + timing)

deseq = DESeq(deseq)
res = results(deseq, contrast = c("condition", "controle", paste("RNAi",condition)))

        png(paste0("./Analyse_petit/", condition,"_",type, "_DESeq_graph",".png"), width = 960, height = 480)
          par(mfrow=c(1,2))
          plotDispEsts(deseq, ylim = c(1e-6, 1e1))
          plotMA(res, ylim = c(-1, 1))
        dev.off()

DEgenes = res[which(as.data.frame(res)$log2FoldChange > FC),]
DEgenes$ID = row.names(DEgenes)

DEanotation = merge(annotation, as.data.frame(DEgenes), by.x = "ID", by.y = "ID")

write.table(DEanotation, paste0("./Analyse_petit/", condition,"_",type, "_DEgenes_DEseq.tab"), sep = "\t", quote = F, row.names = F, col.names = T)


#### Faire un tableau avec les résultats des différentes analyses ####

ctip = read.table(paste0("./Analyse_petit/CTIP_",type, "_DEgenes_DEseq.tab"), header=T,sep="\t")
xrcc4 = read.table(paste0("./Analyse_petit/XRCC4_",type, "_DEgenes_DEseq.tab"), header=T,sep="\t")
both = merge(ctip[,c(1,7,11)], xrcc4[,c(1,7,11)], by.x = "ID", by.y = "ID")
both = both[order(both$SYNONYMS),]
colnames(both)[c(6:9)]= c("CTIP_baseMean","CTIP_log2FoldChange","XRCC4_baseMean", "XRCC4_log2FoldChange")

write.table(both, paste0("./Analyse_petit/CTIPetXRCC4_",type, "_DEgenes_DEseq.tab"), sep = "\t", quote = F, row.names = F, col.names = T)
