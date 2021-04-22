source("2_Mise_en_forme_des_donnees.R")

##### Chargement des librairies pour l'analyse #############
library(FactoMineR)
library(factoextra)
library(DESeq2)

##### Analyse DESeq2 ###########

# Filtration des gènes avec trop peu de compatage (seuil arbitraire)
# countdata = countdata[rowSums(countdata) > 50,]

# Mise en forme des données
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Feeding + Timing)


# Analyse DESeq2
deseq = DESeq(deseq)

###############################################
# Reprise des variables et analyses d'Olivier #
###############################################
notAllZero = (rowSums(counts(deseq)) > 0 )
labels=colnames(countdata)


#### Moyenne des valeurs de comptage normalisées pour chaque point du timining####
##EARLY, INTERMEDIATE, LATE ##
time_points = c(paste("CTRL", timing_ctrl, sep = "_" ),paste(condition, timing_rnai, sep = "_" ))

# Extraction des données de comptage de DESeq2
geneNormCountsTable=counts(deseq,normalized=T)
meanGeneNormCountsTable =data.frame(ID=rownames(geneNormCountsTable))

# Calcule des moyennes ==> Pourquoi utiliser cette partie du code au lieu de la fonction collapseReplicates de DESeq2 ?
for(p in unique(time_points)) {
  
  if(length(labels[time_points==p])==1)  { 
    meanGeneNormCountsTable[,p] =geneNormCountsTable[, labels[time_points==p]]
  } else { 
    meanGeneNormCountsTable[,p] =apply(geneNormCountsTable[, labels[time_points==p]],1,mean) 
  }
  
}
rownames(meanGeneNormCountsTable)=meanGeneNormCountsTable$ID
meanGeneNormCountsTable=meanGeneNormCountsTable[,-1]

#### Comparaison point par point des différents timing ####
comparisons = list(
  "VEG" = time_points[grep("VEG",time_points)],
  "EARLY" = time_points[grep("EARLY",time_points)],
  "INTER" = time_points[grep("INTER",time_points)],
  "LATE" = time_points[grep("LATE",time_points)]
)


significant_up=list()
significant_down=list()

par(mfrow=c(1,4))

## Comparerle controle au RNAi

resContrast=results(deseq,contrast=c("Feeding","controle", paste("RNAi", condition)))
# Affichage résultat 
#summary(resContrast)
#as.data.frame(resContrast) 

resContrast=resContrast[notAllZero,] #Suppression des gènes sans comptages

#Selection des lignes non vide dont la pvalue est inferieur à la pvalue seuil
resContrast_sig = resContrast[!is.na(resContrast$padj) & resContrast$padj < pvalue , ] 
#Selection des lignes dont les log2foldchange sont compris dans le seuil FC appliqué
resContrast_sig = resContrast_sig[resContrast_sig$log2FoldChange >= log2(FC) | resContrast_sig$log2FoldChange <= log2(1/FC), ]

#Extraction du tableau avec les données des gènes considérés dérégulés
resContrast_sig = as.data.frame(resContrast_sig) 

############ Zone de test de FC adapté aux données CTIP #################
# manuel_non_DEG=c("PTET.51.1.G0660118", "PTET.51.1.G1790042", "PTET.51.1.G0030302","PTET.51.1.G1740049","PTET.51.1.G0350134","PTET.51.1.G0230191","PTET.51.1.G0210241","PTET.51.1.G0770102","PTET.51.1.G0220178", "PTET.51.1.G0990073", "PTET.51.1.G0480035", "PTET.51.1.G0370168", "PTET.51.1.G0360062", "PTET.51.1.G1460025", "PTET.51.1.G0170355", "PTET.51.1.G0250220", "PTET.51.1.G0950175", "PTET.51.1.G0170354", "PTET.51.1.G0460033", "PTET.51.1.G1110086", "PTET.51.1.G1280115", "PTET.51.1.G1510135", "PTET.51.1.G1630015", "PTET.51.1.G0640197")
# manuel_DEG_down = c("PTET.51.1.G0540024","PTET.51.1.G0900102","PTET.51.1.G0020380","PTET.51.1.G0120328","PTET.51.1.G0870035","PTET.51.1.G0360066","PTET.51.1.G1010039","PTET.51.1.G0150242","PTET.51.1.G0980137","PTET.51.1.G0170233","PTET.51.1.G0120245","PTET.51.1.G0380022","PTET.51.1.G0230222","PTET.51.1.G0920155","PTET.51.1.G0490162","PTET.51.1.G0070121","PTET.51.1.G0610198","PTET.51.1.G1200062","PTET.51.1.G0680113","PTET.51.1.G1300067","PTET.51.1.G1400105","PTET.51.1.G0030168","PTET.51.1.G0380073","PTET.51.1.G0450225","PTET.51.1.G0620215","PTET.51.1.G0050231","PTET.51.1.G0020217","PTET.51.1.G0350154","PTET.51.1.G0350166","PTET.51.1.G0240239","PTET.51.1.G0010374","PTET.51.1.G0480099","PTET.51.1.G0080368","PTET.51.1.G0210235","PTET.51.1.G0110289","PTET.51.1.G1530110","PTET.51.1.G0260051","PTET.51.1.G1150114","PTET.51.1.G0590028","PTET.51.1.G0340197","PTET.51.1.G0110267","PTET.51.1.G1140146","PTET.51.1.G0350167","PTET.51.1.G0370136","PTET.51.1.G0210213","PTET.51.1.G0220140","PTET.51.1.G0360089","PTET.51.1.G0010451","PTET.51.1.G0020335","PTET.51.1.G0250013","PTET.51.1.G1330044","PTET.51.1.G0060034","PTET.51.1.G0530071")
# manuel_non_DEG = as.data.frame(resContrast)[is.element(row.names(as.data.frame(resContrast)),manuel_non_DEG),]
# summary(manuel_non_DEG)
# manuel_DEG_down = as.data.frame(resContrast)[is.element(row.names(as.data.frame(resContrast)), manuel_DEG_down),]
# summary(manuel_DEG_down)
# 
# log2FC = list(manuel_non_DEG$log2FoldChange, manuel_DEG_down$log2FoldChange)
# names(log2FC) =c("Non DEG", "DEG")
# boxplot(log2FC, main = paste("log2FC","\n",i), horizontal=F)
##############