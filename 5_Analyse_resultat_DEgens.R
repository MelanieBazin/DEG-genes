options(stringsAsFactors = FALSE)
library("stringr") 
library(ggvenn)
library(ggplot2)
library("RColorBrewer")

source("0_Cluster.R")

# Definitir les fichiers à ouvrir
date = "02-08"
condition =  names(rnai_list)[2]

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
path = paste0("./Analyse/",file_name, "/", condition, "/")
save_path = paste0(path, "Analyse/")
dir.create(save_path,recursive=T,showWarnings=F)

# Donnée externe
TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation = annotation[,c(1,3:5,13,6:11,2)]
rownames(annotation)=annotation$ID

# Donnée DESeq2
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]

#### Ouverture tableau de données des gènes DEG ####
TAB = list()
for (R in RNAi){
  if (R == "CTIP"){
    TAB = c(TAB, CTIP_early = list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_EARLY_NoFilter.tab"), header=T,sep="\t",quote='')))
    TAB = c(TAB, CTIP_inter = list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_INTER_NoFilter.tab"), header=T,sep="\t",quote='')))
  }else {
    TAB = c(TAB, list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_LATE_NoFilter.tab"), header=T,sep="\t",quote='')))
  }
}
names(TAB)[3:length(names(TAB))]= RNAi[-1]

# Definir les sous-liste de genes
source("5-1_Filtres.R")

#### Summary table ####
summary_tab = annotation
for (cond in names(TAB)){
  tab = TAB[[cond]]
  mini_tab = cbind(tab$ID, tab$log2FoldChange, tab$padj, tab$REGULATION)
  colnames(mini_tab) = c("ID", paste(cond, c("log2FC", "padj", "REG"), sep = "_"))
  summary_tab = merge(summary_tab, mini_tab, by = "ID", all = T)
}
colnames(TurboPGM) = c("PROTEIN_NAME", paste0("TurboPGM_", c("log2FC", "-log10pval")))
colnames(TurboPGML4) = c("PROTEIN_NAME", paste0("TurboPGML4_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGM, by = "PROTEIN_NAME", all = T)

write.table(summary_tab,paste0("Analyse/",file_name,"/",condition,"/Summary_",condition,".tab"), sep = "\t", row.names = F) 

### Venn Diagrame ####
###### Croisement des tubo ####
LIST = list(
  TurboPGM = TurboPGM$PROTEIN_NAME,
  TurboPGML4 = TurboPGML4$PROTEIN_NAME
)
png(paste0(save_path,"Venn_Turbo.png"))
ggvenn(LIST,
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7)
dev.off()

turboPGM = annotation$ID[which(is.element(annotation$PROTEIN_NAME, TurboPGM$PROTEIN_NAME))]
turboPGML4 = annotation$ID[which(is.element(annotation$PROTEIN_NAME, TurboPGML4$PROTEIN_NAME))]
turbo = intersect(turboPGM, turboPGML4)

###### Croisement des UP PGM KU80c XRCC4 ####
LIST = UP_PKX
png(paste0(save_path,"Venn_UP_PKX.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

up_pkx = intersect(UP_PKX[["UP_XRCC4"]], intersect(UP_PKX[["UP_PGM"]],UP_PKX[["UP_KU80c"]]))

# Avec intermediate peak
LIST = c(UP_PKX, Intermediate_peak = list(inter_genes))
png(paste0(save_path,"Venn_UP_PKX_INTER.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec intermediate et turbo
LIST = list(UP_PKX = up_pkx,
            TurboPGM = turboPGM,
            TurboPGML4 = turboPGML4,
            Intermediate_peak = inter_genes)
png(paste0(save_path,"Venn_PKX_INTER_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec early peak
LIST = c(UP_PKX, Early_peak = list(early_genes))
png(paste0(save_path,"Venn_UP_PKX_EARLY.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec early et turbo
LIST = list(UP_PKX = up_pkx,
            TurboPGM = turboPGM,
            TurboPGML4 = turboPGML4,
            Early_peak = early_genes)
png(paste0(save_path,"Venn_PKX_EARLY_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

###### Croisement des DOWN CTIP ####
# Avec les données TurboID
LIST = list(DOWN_CTIP = DOWN_C,
            TurboPGM = turboPGM,
            TurboPGML4 = turboPGML4)
png(paste0(save_path,"Venn_CTIP_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec les UP dérégulés en PGM, KU80c, XRCC4
LIST = c(UP_PKX, DOWN_CTIP = list(DOWN_C))
png(paste0(save_path,"Venn_CTIP_PKX.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec les UP dérégulés + intermediate peak
LIST = list(UP_PKX = up_pkx,
            DOWN_CTIP = DOWN_C,
            Intermediate_peak = inter_genes)
png(paste0(save_path,"Venn_CTIP_PKX_INTER.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec les UP dérégulés + intermediate peak + turbo
LIST = list(UP_PKX = up_pkx,
            DOWN_CTIP = DOWN_C,
            TurboPGM_PGML4 = turbo,
            Intermediate_peak = inter_genes)
png(paste0(save_path,"Venn_CTIP_PKX_INTER_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec les UP dérégulés + early peak
LIST = list(UP_PKX = up_pkx,
            DOWN_CTIP = DOWN_C,
            Intermediate_peak = inter_genes)
png(paste0(save_path,"Venn_CTIP_PKX_EARLY.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec les UP dérégulés + early peak + Turbo
LIST = list(UP_PKX = up_pkx,
            DOWN_CTIP = DOWN_C,
            TurboPGM_PGML4 = turbo,
            Intermediate_peak = inter_genes)
png(paste0(save_path,"Venn_CTIP_PKX_EARLY_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()


### Barplot classes des gènes ####
row_order = c("Early peak", "Intermediate peak", "Late peak", "Early repression" ,"Late induction", "Late repression", "none" )
colors = c("purple3","red2","chartreuse4","dodgerblue3","deeppink","darkorange","snow3")

###### Sur UP PGM KU80c & XRCC4 ####
profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
for (n in names(UP_PKX)){
  tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, UP_PKX[[n]]))]))
  profil = merge(profil, tab, by = "Var1", all = T)
  
}
tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, up_pkx))]))
profil = merge(profil, tab, by = "Var1", all = T)

rownames(profil) = profil$Var1
profil = profil[,-1]
colnames(profil) = c("ALL", names(UP_PKX), "UP_ALL")

# Réordonner les lignes
profil = profil[row_order,]


### Histogramme empilés
png(paste0(save_path,"Profils_barplot_UP.png"),width = 550, height = 500)
barplot(as.matrix(profil),
        col = colors,
        main = "Profil repartition of UP deregulated genes",
        ylab = "gene nb")

legend("topright",
       legend = rownames(profil),
       fill = colors,
       bty = "n")
dev.off()

# Création d'un tableau avec ses pourcentages
profil_prct = profil
for (n in 1:ncol(profil)){
  profil_prct[,n] = profil_prct[,n]/sum(profil[,n])*100
}

png(paste0(save_path,"Profils_barplot_UP_prct.png"),width = 550, height = 500)
barplot(as.matrix(profil_prct),
        col = colors,
        main = "Profil repartition of UP deregulated genes",
        ylab = "% of genes",
        names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
dev.off()

### Histogramme enrichissement

###### Sur DOWN CTIP + UP PKX ####
# Histogramme empilés

# Histogramme enrichissement

### Histogrammes IES in genes ####
###### Sur UP PGM KU80c & XRCC4 ####
# Histogramme empilés

# Histogramme enrichissement

###### Sur DOWN CTIP + UP PKX ####
# Histogramme empilés

# Histogramme enrichissement




