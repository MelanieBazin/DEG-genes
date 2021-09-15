# install.packages("RColorBrewer")
# install.packages("pheatmap")
# install.packages("seqinr")
# install.packages("latticeExtra") # Erreur d'installation sur R4.0.4
# install.packages("FactoMineR")
# install.packages("factoextra")
# install.packages("gplots")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("viridis")
# install.packages("caret")
# install.packages("magick")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("vsn")
# BiocManager::install("DESeq")
# BiocManager::install("DESeq2")
# BiocManager::install("topGO")
# BiocManager::install("latticeExtra") # Erreur d'installation sur R4.0.4
# BiocManager::install("lifecycle")
# BiocManager::install("heatmaps") # Erreur d'installation sur R4.0.4
# BiocManager::install("sva")
# BiocManager::install("klaR")
# BiocManager::install("ComplexHeatmap")
# 
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")



options(stringsAsFactors = FALSE)
library("stringr") 

annotation = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
annotation$NAME[grep("PTET.51.1.G0110267",annotation$ID)]="PGML1"
annotation$NAME[grep("PTET.51.1.G0380073",annotation$ID)]="PGML2"
annotation$NAME[grep("PTET.51.1.G0010374",annotation$ID)]="PGML3a"
annotation$NAME[grep("PTET.51.1.G0080308",annotation$ID)]="PGML3b"
annotation$NAME[grep("PTET.51.1.G0020217",annotation$ID)]="PGML3c"
annotation$NAME[grep("PTET.51.1.G0340197",annotation$ID)]="PGML4a"
annotation$NAME[grep("PTET.51.1.G0480099",annotation$ID)]="PGML4b"
annotation$NAME[grep("PTET.51.1.G0570051",annotation$ID)]="PGML5a"
annotation$NAME[grep("PTET.51.1.G0510172",annotation$ID)]="PGML5b"

annotation$NAME[grep("PTET.51.1.G1110086",annotation$ID)]="XRCC4"
annotation$NAME[grep("PTET.51.1.G0020380",annotation$ID)]="CERa"
annotation$NAME[grep("PTET.51.1.G0220178",annotation$ID)]="CERb"
annotation$NAME[grep("PTET.51.1.G0540024",annotation$ID)]="LigIVa"
annotation$NAME[grep("PTET.51.1.G0720019",annotation$ID)]="LigIVb"

annotation$NAME[grep("PTET.51.1.G1330044",annotation$ID)]="EZL2"
annotation$NAME[grep("PTET.51.1.G1740049",annotation$ID)]="EZL1"

annotation$NAME[grep("PTET.51.1.G0490126",annotation$ID)]=" "
annotation$NAME[grep("PTET.51.1.G0640198",annotation$ID)]="Mre11a?"

annotation$NAME[grep("PTET.51.1.G0380048",annotation$ID)]="SPT5v"
annotation$NAME[grep("PTET.51.1.G0770102",annotation$ID)]="SPT5m"

annotation$NAME[grep("PTET.51.1.G0780031",annotation$ID)]="CAF1"
annotation$NAME[grep("PTET.51.1.G0230191",annotation$ID)]="Spo11"

annotation$NAME[grep("PTET.51.1.G0210235",annotation$ID)]="PolXa"
annotation$NAME[grep("PTET.51.1.G0360066",annotation$ID)]="PolXb"
annotation$NAME[grep("PTET.51.1.G0460033",annotation$ID)]="PolXc"
annotation$NAME[grep("PTET.51.1.G1010039",annotation$ID)]="PolXd"


annotation$NAME[grep("PTET.51.1.G0330337",annotation$ID)]="HMG2a"
annotation$NAME[grep("PTET.51.1.G0450312",annotation$ID)]="HMG2b"
annotation$NAME[grep("PTET.51.1.G0380231",annotation$ID)]="HMG2c"
annotation$NAME[grep("PTET.51.1.G0500224",annotation$ID)]="HMG2d"

annotation$NAME[grep("PTET.51.1.G0550242",annotation$ID)]="HMG3a"
annotation$NAME[grep("PTET.51.1.G1080160",annotation$ID)]="HMG3b"
annotation$NAME[grep("PTET.51.1.G1210130",annotation$ID)]="HMG3c"
annotation$NAME[grep("PTET.51.1.G1050131",annotation$ID)]="HMG3d"
annotation$NAME[grep("PTET.51.1.G1450030",annotation$ID)]="HMG3e"
annotation$NAME[grep("PTET.51.1.G1520034",annotation$ID)]="HMG3f"
annotation$NAME[grep("PTET.51.1.G0800204",annotation$ID)]="HMG3g"

annotation$NAME[grep("PTET.51.1.G1280139",annotation$ID)]="HMG6"

annotation$NAME[grep("PTET.51.1.G1660041",annotation$ID)]="HMG7a"
annotation$NAME[grep("PTET.51.1.G0400270",annotation$ID)]="HMG7b"
annotation$NAME[grep("PTET.51.1.G1100130",annotation$ID)]="HMG7c"
annotation$NAME[grep("PTET.51.1.G1110064",annotation$ID)]="HMG7e"
annotation$NAME[grep("PTET.51.1.G0850121",annotation$ID)]="HMG7f"
annotation$NAME[grep("PTET.51.1.G1420074",annotation$ID)]="HMG7g"

annotation$NAME[grep("PTET.51.1.G1550114",annotation$ID)]="HMG10a"
annotation$NAME[grep("PTET.51.1.G1070037",annotation$ID)]="HMG10b"

annotation$NAME[grep("PTET.51.1.G0540171",annotation$ID)]="HMG13a"
annotation$NAME[grep("PTET.51.1.G1350008",annotation$ID)]="HMG13b"
annotation$NAME[grep("PTET.51.1.G0390142",annotation$ID)]="HMG13d"

annotation$NAME[grep("PTET.51.1.G0100092",annotation$ID)]="HMG14a"
annotation$NAME[grep("PTET.51.1.G0050352",annotation$ID)]="HMG14c"

annotation$NAME[grep("PTET.51.1.G1490060",annotation$ID)]="HMG15a"
annotation$NAME[grep("PTET.51.1.G1270097",annotation$ID)]="HMG15b"
annotation$NAME[grep("PTET.51.1.G0560078",annotation$ID)]="HMG15c"
annotation$NAME[grep("PTET.51.1.G0580048",annotation$ID)]="HMG15d"

annotation$NAME[grep("PTET.51.1.G0090176",annotation$ID)]="HMG18a"
annotation$NAME[grep("PTET.51.1.G0260144",annotation$ID)]="HMG18b"
annotation$NAME[grep("PTET.51.1.G0650159",annotation$ID)]="HMG18c"


mt = read.table("./DATA/mtF_ID.csv", sep = ";", header = T)
for (i in mt$ID){
  annotation$NAME[grep(i,annotation$ID)]=mt$SYNONYMS[grep(i,mt$ID)]
}

for (i in grep("PTIWI",annotation$SYNONYMS)){
  annotation$NAME[i]= str_split(annotation$SYNONYMS[i],",")[[1]][1]
}

# gene_autogamy=read.table("DATA/autogamy_ptetraurelia_mac_51_annotation_v2.0_significant.tab", header=T, sep="\t")
# annotation = merge(annotation[,1:3],gene_autogamy[,c(1,6)],by.x = "ID", by.y = "ID")
# deg_eva = read.table("DATA/EVA_Siginificantly_DEGenes_PGM_KU80.tab", header=T, sep="\t")
# annotation$PGM = NA
# annotation$KU80 = NA
# for (z in 1:nrow(deg_eva)){
#   annotation$PGM[grep(deg_eva$ID[z],annotation$ID)]=deg_eva$PGM_REGULATION[z]
#   annotation$KU80[grep(deg_eva$ID[z],annotation$ID)]=deg_eva$KU80_REGULATION[z]
# }


# gene_ies=read.delim("../DATA/PARAMECIUM/GENOMIC/tetraurelia/macronucleus/ANNOTATION/ptetraurelia_mac_51/v2/ptetraurelia_mac_51_annotation_v2.0.ies",h=T)
# annotation = merge(gene_ies[,1:2],annotation, by = "ID")
annotation1 = read.table("./DATA/My_annotation.tab", header = T, sep = "\t" )
annotation = merge(annotation1[,1:2],annotation, by = "ID")



write.table(annotation, "./DATA/My_annotation2.tab", sep = "\t", row.names = F)


