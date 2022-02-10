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

# Definir les sous-liste de genes
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]

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
path = paste0(save_path,"VennDiagrame/")
dir.create(path ,recursive=T,showWarnings=F)
###### Croisement des tubo ####
LIST = list(
  TurboPGM = TurboPGM$PROTEIN_NAME,
  TurboPGML4 = TurboPGML4$PROTEIN_NAME
)
png(paste0(path,"Venn_Turbo.png"))
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
png(paste0(path,"Venn_UP_PKX.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec intermediate peak
LIST = c(UP_PKX, Intermediate_peak = list(inter_genes))
png(paste0(path,"Venn_UP_PKX_INTER.png"))
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
png(paste0(path,"Venn_PKX_INTER_TURBO.png"))
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
png(paste0(path,"Venn_UP_PKX_EARLY.png"))
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
png(paste0(path,"Venn_PKX_EARLY_TURBO.png"))
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
png(paste0(path,"Venn_CTIP_TURBO.png"))
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
png(paste0(path,"Venn_CTIP_PKX.png"))
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
png(paste0(path,"Venn_CTIP_PKX_INTER.png"))
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
png(paste0(path,"Venn_CTIP_PKX_INTER_TURBO.png"))
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
png(paste0(path,"Venn_CTIP_PKX_EARLY.png"))
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
png(paste0(path,"Venn_CTIP_PKX_EARLY_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()


### Barplot classes des gènes ####
path = paste0(save_path,"Barplot_profil/")
dir.create(path ,recursive=T,showWarnings=F)

row_order = c("Early peak", "Intermediate peak", "Late peak", "Early repression" ,"Late induction", "Late repression", "none" )
colors = c("purple3","red2","chartreuse4","dodgerblue3","deeppink","darkorange","snow3")

###### Sur UP PGM KU80c & XRCC4 ####
UP_PKX = c(UP_PKX, UP_ALL = list(up_pkx))

profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
for (n in names(UP_PKX)){
  tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, UP_PKX[[n]]))]))
  profil = merge(profil, tab, by = "Var1", all = T)
  
}

rownames(profil) = profil$Var1
profil = profil[,-1]
colnames(profil) = c("ALL", names(UP_PKX))

# Réordonner les lignes
profil = profil[row_order,]


### Histogramme empilés
png(paste0(path,"Profils_barplot_UP.png"),width = 550, height = 500)
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

png(paste0(path,"Profils_barplot_UP_prct.png"),width = 550, height = 500)
barplot(as.matrix(profil_prct),
        col = colors,
        main = "Profil repartition of UP deregulated genes",
        ylab = "% of genes",
        names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
dev.off()

### Histogramme enrichissement
profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
for(up in names(UP_PKX)){
  tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, UP_PKX[[up]]))]))
  tab = merge(profil, tab, by = "Var1")
  rownames(tab) = tab[,"Var1"]
  tab = tab[,-1]
  colnames(tab) = c("ALL",up)
  
  tab_prct = tab
  tab_prct[,"ALL"] = tab[,"ALL"]/length(annotation$ID)*100
  for (n in rownames(tab)){
    tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
  }
  
  
  png(paste0(path,"Profils_barplot_",up,".png"),width = 800, height = 500)
  barplot(t(as.matrix(tab_prct)),
          beside = T,
          main = "Profil repartition of UP deregulated genes",
          ylab = "% of genes",
          ylim = c(0,60),
          col = c("grey", "indianred2"),
          names.arg = sub(" "," \n ",rownames(tab)))
  
  legend("topleft",
         legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
         fill = c("grey", "indianred2"),
         bty = "n")
  dev.off()
}

###### Sur DOWN CTIP + UP PKX ####
stdCTIP = list(DOWN_CTIP = DOWN_C,
            UP_ALL = up_pkx,
            DOWN_UP = intersect(DOWN_C, up_pkx))

profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
for (n in names(stdCTIP)){
  tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, stdCTIP[[n]]))]))
  profil = merge(profil, tab, by = "Var1", all = T)
  
}

rownames(profil) = profil$Var1
profil = profil[,-1]
colnames(profil) = c("ALL", names(stdCTIP))

# Réordonner les lignes
profil = profil[row_order,]


### Histogramme empilés
png(paste0(path,"Profils_barplot_CTIP.png"),width = 450, height = 500)
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

png(paste0(path,"Profils_barplot_CTIP_prct.png"),width = 450, height = 500)
barplot(as.matrix(profil_prct),
        col = colors,
        main = "Profil repartition of UP deregulated genes",
        ylab = "% of genes",
        names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
dev.off()

### Histogramme enrichissement
profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
for(up in names(stdCTIP)){
  tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, stdCTIP[[up]]))]))
  tab = merge(profil, tab, by = "Var1")
  rownames(tab) = tab[,"Var1"]
  tab = tab[,-1]
  colnames(tab) = c("ALL",up)
  
  tab_prct = tab
  tab_prct[,"ALL"] = tab[,"ALL"]/length(annotation$ID)*100
  for (n in rownames(tab)){
    tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
  }
  
  png(paste0(path,"Profils_barplot_",up,".png"),width = 800, height = 500)
  barplot(t(as.matrix(tab_prct)),
          beside = T,
          main = "Profil repartition of UP deregulated genes",
          ylab = "% of genes",
          ylim = c(0,60),
          col = c("grey", "indianred2"),
          names.arg = sub(" "," \n ",rownames(tab)))
  
  legend("topleft",
         legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
         fill = c("grey", "indianred2"),
         bty = "n")
  dev.off()
}

### Histogrammes IES in genes ####
path = paste0(save_path,"Barplot_IES/")
dir.create(path ,recursive=T,showWarnings=F)

colors = c("darkblue","snow3")

###### Sur UP PGM KU80c & XRCC4 ####
profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                         sum(annotation$NB_IES == 0)), 
                       row.names = c("IES+", "IES-"))
for (n in names(UP_PKX)){
  tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, UP_PKX[[n]]))] != 0),
          sum(annotation$NB_IES[which(is.element(annotation$ID, UP_PKX[[n]]))] == 0))
  profil = cbind(profil, tab)
  
}

colnames(profil) = c("ALL", names(UP_PKX))


### Histogramme empilés
png(paste0(path,"Profils_barplot_UP.png"),width = 550, height = 500)
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

png(paste0(path,"Profils_barplot_UP_prct.png"),width = 550, height = 500)
barplot(as.matrix(profil_prct),
        col = colors,
        main = "Profil repartition of UP deregulated genes",
        ylab = "% of genes",
        names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
dev.off()

### Histogramme enrichissement
profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                         sum(annotation$NB_IES == 0)), 
                       row.names = c("IES+", "IES-"))
for(up in names(UP_PKX)){
  tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, UP_PKX[[up]]))] != 0),
          sum(annotation$NB_IES[which(is.element(annotation$ID, UP_PKX[[up]]))] == 0))
  tab = cbind(profil, tab)
  colnames(tab) = c("ALL",up)
  
  tab_prct = tab
  tab_prct[,"ALL"] = tab[,"ALL"]/length(annotation$ID)*100
  for (n in rownames(tab)){
    tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
  }
  
  
  png(paste0(path,"Profils_barplot_",up,".png"),width = 400, height = 500)
  barplot(t(as.matrix(tab_prct)),
          beside = T,
          main = "Profil repartition of UP deregulated genes",
          ylab = "% of genes",
          ylim = c(0,60),
          col = c("grey", "indianred2"),
          names.arg = sub(" "," \n ",rownames(tab)))
  
  legend("topleft",
         legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
         fill = c("grey", "indianred2"),
         bty = "n")
  dev.off()
}

###### Sur DOWN CTIP + UP PKX ####
profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                         sum(annotation$NB_IES == 0)), 
                       row.names = c("IES+", "IES-"))
for (n in names(stdCTIP)){
  tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, stdCTIP[[n]]))] != 0),
          sum(annotation$NB_IES[which(is.element(annotation$ID, stdCTIP[[n]]))] == 0))
  profil = cbind(profil, tab)
  
}

colnames(profil) = c("ALL", names(stdCTIP))


### Histogramme empilés
png(paste0(path,"Profils_barplot_UP.png"),width = 550, height = 500)
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

png(paste0(path,"Profils_barplot_UP_prct.png"),width = 550, height = 500)
barplot(as.matrix(profil_prct),
        col = colors,
        main = "Profil repartition of UP deregulated genes",
        ylab = "% of genes",
        names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
dev.off()

### Histogramme enrichissement
profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                         sum(annotation$NB_IES == 0)), 
                       row.names = c("IES+", "IES-"))
for(up in names(stdCTIP)){
  tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, stdCTIP[[up]]))] != 0),
          sum(annotation$NB_IES[which(is.element(annotation$ID, stdCTIP[[up]]))] == 0))
  tab = cbind(profil, tab)
  colnames(tab) = c("ALL",up)
  
  tab_prct = tab
  tab_prct[,"ALL"] = tab[,"ALL"]/length(annotation$ID)*100
  for (n in rownames(tab)){
    tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
  }
  
  
  png(paste0(path,"Profils_barplot_",up,".png"),width = 400, height = 500)
  barplot(t(as.matrix(tab_prct)),
          beside = T,
          main = "Profil repartition of UP deregulated genes",
          ylab = "% of genes",
          ylim = c(0,60),
          col = c("grey", "indianred2"),
          names.arg = sub(" "," \n ",rownames(tab)))
  
  legend("topleft",
         legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
         fill = c("grey", "indianred2"),
         bty = "n")
  dev.off()
}

sink(paste0(save_path,"/sessionInfo.txt"))
print(sessionInfo())
sink()


