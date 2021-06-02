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
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("vsn")
# BiocManager::install("DESeq")
# BiocManager::install("DESeq2")
# BiocManager::install("topGO")
# BiocManager::install("latticeExtra")
# BiocManager::install("lifecycle")
# BiocManager::install("heatmaps")
# BiocManager::install("sva")
# BiocManager::install("klaR")
# 
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")



options(stringsAsFactors = FALSE)

annotation = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
annotation$NAME[grep("PTET.51.1.G0490162",annotation$ID)]="PGM, PiggyMac"
annotation$NAME[grep("PTET.51.1.G0110267",annotation$ID)]="PGML1"
annotation$NAME[grep("PTET.51.1.G0380073",annotation$ID)]="PGML2"
annotation$NAME[grep("PTET.51.1.G0010374",annotation$ID)]="PGML3a"
annotation$NAME[grep("PTET.51.1.G0080308",annotation$ID)]="PGML3b"
annotation$NAME[grep("PTET.51.1.G0020217",annotation$ID)]="PGML3c"
annotation$NAME[grep("PTET.51.1.G0340197",annotation$ID)]="PGML4a"
annotation$NAME[grep("PTET.51.1.G0480099",annotation$ID)]="PGML4b"
annotation$NAME[grep("PTET.51.1.G0570051",annotation$ID)]="PGML5a"
annotation$NAME[grep("PTET.51.1.G0510172",annotation$ID)]="PGML5b"

annotation$NAME[grep("PTET.51.1.G0150242",annotation$ID)]="Ku70a"
annotation$NAME[grep("PTET.51.1.G0250220",annotation$ID)]="Ku70b"
annotation$NAME[grep("PTET.51.1.G1460025",annotation$ID)]="Ku80a"
annotation$NAME[grep("PTET.51.1.G1510135",annotation$ID)]="Ku80b"
annotation$NAME[grep("PTET.51.1.G1140146",annotation$ID)]="Ku80c"

annotation$NAME[grep("PTET.51.1.G1110086",annotation$ID)]="XRCC4"
annotation$NAME[grep("PTET.51.1.G0020380",annotation$ID)]="CERa"
annotation$NAME[grep("PTET.51.1.G0220178",annotation$ID)]="CERb"
annotation$NAME[grep("PTET.51.1.G0540024",annotation$ID)]="LigIVa"
annotation$NAME[grep("PTET.51.1.G0720019",annotation$ID)]="LigIVb"

annotation$NAME[grep("PTET.51.1.G1330044",annotation$ID)]="EZL2"
annotation$NAME[grep("PTET.51.1.G1740049",annotation$ID)]="EZL1"

annotation$NAME[grep("PTET.51.1.G0490126",annotation$ID)]=" "
annotation$NAME[grep("PTET.51.1.G0640197",annotation$ID)]="Mre11a"
annotation$NAME[grep("PTET.51.1.G0640198",annotation$ID)]="Mre11a?"
annotation$NAME[grep("PTET.51.1.G0790183",annotation$ID)]="Mre11b"

annotation$NAME[grep("PTET.51.1.G0380048",annotation$ID)]="SPT5v"
annotation$NAME[grep("PTET.51.1.G0770102",annotation$ID)]="SPT5m"

annotation$NAME[grep("PTET.51.1.G0650078",annotation$ID)]="CtIPa"
annotation$NAME[grep("PTET.51.1.G0980137",annotation$ID)]="CtIPb"

annotation$NAME[grep("PTET.51.1.G0780031",annotation$ID)]="CAF1"
annotation$NAME[grep("PTET.51.1.G0230191",annotation$ID)]="Spo11"

annotation$NAME[grep("PTET.51.1.G0210235",annotation$ID)]="PolXa"
annotation$NAME[grep("PTET.51.1.G0360066",annotation$ID)]="PolXb"
annotation$NAME[grep("PTET.51.1.G0460033",annotation$ID)]="PolXc"
annotation$NAME[grep("PTET.51.1.G1010039",annotation$ID)]="PolXd"

mtF = read.csv("DATA/mtF_ID.csv", header=T, sep=";")
for (i in 1:nrow(mtF)){
  annotation$SYNONYMS[grep(mtF$ID[i],annotation$ID)]=mtF$SYNONYMS[i]
}

gene_autogamy=read.table("DATA/autogamy_ptetraurelia_mac_51_annotation_v2.0_significant.tab", header=T, sep="\t")
annotation = merge(annotation[,1:3],gene_autogamy[,c(1,6)],by.x = "ID", by.y = "ID")
deg_eva = read.table("DATA/EVA_Siginificantly_DEGenes_PGM_KU80.tab", header=T, sep="\t")
annotation$PGM = NA
annotation$KU80 = NA
for (z in 1:nrow(deg_eva)){
  annotation$PGM[grep(deg_eva$ID[z],annotation$ID)]=deg_eva$PGM_REGULATION[z]
  annotation$KU80[grep(deg_eva$ID[z],annotation$ID)]=deg_eva$KU80_REGULATION[z]
}


# gene_ies=read.delim("../DATA/PARAMECIUM/GENOMIC/tetraurelia/macronucleus/ANNOTATION/ptetraurelia_mac_51/v2/ptetraurelia_mac_51_annotation_v2.0.ies",h=T)
# annotation = merge(gene_ies[,1:2],annotation, by = "ID")
annotation1 = read.table("./DATA/My_annotation.tab", header = T, sep = "\t" )
annotation = merge(annotation1[,1:2],annotation, by = "ID")



write.table(annotation, "./DATA/My_annotation2.tab", sep = "\t", row.names = F)


