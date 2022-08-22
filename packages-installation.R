# install.packages("RColorBrewer")
# install.packages("pheatmap")
# install.packages("seqinr")
# install.packages("latticeExtra")
# install.packages("FactoMineR")
# install.packages("factoextra")
# install.packages("gplots")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("viridis")
# install.packages("caret")
# install.packages("magick")
# install.packages("devtools")
# install.packages("ggseqlogo")
# install.packages("stringr")
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocVersion")
# BiocManager::install("vsn")
# BiocManager::install("DESeq2")
# BiocManager::install("topGO")
# BiocManager::install("latticeExtra") 
# BiocManager::install("lifecycle")
# BiocManager::install("heatmaps") 
# BiocManager::install("sva")
# BiocManager::install("klaR")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("memes")
# 
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
# devtools::install_github("zhangyuqing/sva-devel")
# library(devtools)
# install_github("vqv/ggbiplot")

options(stringsAsFactors = FALSE)
library(stringr) 

annotation = read.table("./DATA/Annotation_ptetraurelia_mac_51_Sherlock20220819.tab",header=T,sep="\t",quote='')

annotation$ID = rownames(annotation)
colnames(annotation)[grep("Number.of.IESs.within.gene",colnames(annotation))] = "NB_IES"
colnames(annotation)[grep("Autogamy.time.course.expression.group",colnames(annotation))] = "EXPRESSION_PROFIL"
annotation = annotation[,c(12,1:4,11,5:10)]

write.table(annotation, "./DATA/My_annotation.tab", sep = "\t", row.names = F)

# Add private annotation

annotation$Name[grep("PTET.51.1.G0720019",annotation$ID)]="LIG4b"

annotation$Name[grep("PTET.51.1.G0640198",annotation$ID)]="Mre11a?"

annotation$Name[grep("PTET.51.1.G0210235",annotation$ID)]="PolXa"
annotation$Name[grep("PTET.51.1.G0360066",annotation$ID)]="PolXb"
annotation$Name[grep("PTET.51.1.G0460033",annotation$ID)]="PolXc"
annotation$Name[grep("PTET.51.1.G1010039",annotation$ID)]="PolXd"


annotation$Name[grep("PTET.51.1.G0330337",annotation$ID)]="HMGB2a"
annotation$Name[grep("PTET.51.1.G0450312",annotation$ID)]="HMGB2b"
annotation$Name[grep("PTET.51.1.G0380231",annotation$ID)]="HMGB2c"
annotation$Name[grep("PTET.51.1.G0500224",annotation$ID)]="HMGB2d"

annotation$Name[grep("PTET.51.1.G0550242",annotation$ID)]="HMGB3a"
annotation$Name[grep("PTET.51.1.G1080160",annotation$ID)]="HMGB3b"
annotation$Name[grep("PTET.51.1.G1210130",annotation$ID)]="HMGB3c"
annotation$Name[grep("PTET.51.1.G1050131",annotation$ID)]="HMGB3d"
annotation$Name[grep("PTET.51.1.G1450030",annotation$ID)]="HMGB3e"
annotation$Name[grep("PTET.51.1.G1520034",annotation$ID)]="HMGB3f"
annotation$Name[grep("PTET.51.1.G0800204",annotation$ID)]="HMGB3g"

annotation$Name[grep("PTET.51.1.G1280139",annotation$ID)]="HMGB6"

annotation$Name[grep("PTET.51.1.G1660041",annotation$ID)]="HMGB7a"
annotation$Name[grep("PTET.51.1.G0400270",annotation$ID)]="HMGB7b"
annotation$Name[grep("PTET.51.1.G1100130",annotation$ID)]="HMGB7c"
annotation$Name[grep("PTET.51.1.G1110064",annotation$ID)]="HMGB7e"
annotation$Name[grep("PTET.51.1.G0850121",annotation$ID)]="HMGB7f"
annotation$Name[grep("PTET.51.1.G1420074",annotation$ID)]="HMGB7g"

annotation$Name[grep("PTET.51.1.G1550114",annotation$ID)]="HMGB10a"
annotation$Name[grep("PTET.51.1.G1070037",annotation$ID)]="HMGB10b"

annotation$Name[grep("PTET.51.1.G0540171",annotation$ID)]="HMGB13a"
annotation$Name[grep("PTET.51.1.G1350008",annotation$ID)]="HMGB13b"
annotation$Name[grep("PTET.51.1.G0390142",annotation$ID)]="HMG13d"

annotation$Name[grep("PTET.51.1.G0100092",annotation$ID)]="HMGB14a"
annotation$Name[grep("PTET.51.1.G0050352",annotation$ID)]="HMGB14c"

annotation$Name[grep("PTET.51.1.G1490060",annotation$ID)]="HMGB15a"
annotation$Name[grep("PTET.51.1.G1270097",annotation$ID)]="HMGB15b"
annotation$Name[grep("PTET.51.1.G0560078",annotation$ID)]="HMGB15c"
annotation$Name[grep("PTET.51.1.G0580048",annotation$ID)]="HMGB15d"

annotation$Name[grep("PTET.51.1.G0090176",annotation$ID)]="HMGB18a"
annotation$Name[grep("PTET.51.1.G0260144",annotation$ID)]="HMGB18b"
annotation$Name[grep("PTET.51.1.G0650159",annotation$ID)]="HMGB18c"

annotation$Name[grep("PTET.51.1.G0620188",annotation$ID)]="PHD-finger_Swart"

annotation$Name[grep("PTET.51.1.G0140243",annotation$ID)]="ISWI1a"
annotation$Name[grep("PTET.51.1.G0420126",annotation$ID)]="ISWI1b"
annotation$Name[grep("PTET.51.1.G0290211",annotation$ID)]="ISWI1c"
annotation$Name[grep("PTET.51.1.G0410183",annotation$ID)]="ISWI1d"

annotation$Name[grep("PTET.51.1.G1470016",annotation$ID)]="ISWI2"

mt = read.table("./DATA/mtF_ID.csv", sep = ";", header = T)
for (i in mt$ID){
  annotation$Name[grep(i,annotation$ID)]=mt$SYNONYMS[grep(i,mt$ID)]
}

for (i in grep("PTIWI",annotation$Aliases)){
  annotation$Name[i]= str_split(annotation$Aliases[i],",")[[1]][1]
}


write.table(annotation, "./DATA/My_annotation2.tab", sep = "\t", row.names = F)


