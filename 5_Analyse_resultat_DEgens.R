options(stringsAsFactors = FALSE)
library(stringr) 
library(ggvenn)
library(ggplot2)
library(RColorBrewer)

source("0_Cluster.R")
source("0_Visualisation_fonction.R")

# Definitir les fichiers à ouvrir
date = Sys.Date()
# date = "2022-02-21"
condition =  names(rnai_list)[2]

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
path = paste0("./Analyse/",file_name, "/", condition, "/")
save_path = paste0(path, "Analyse/")
dir.create(save_path,recursive=T,showWarnings=F)

# Definir les sous-liste de genes
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]

source("5-1_Filtres.R")

#### Summary table ####
print("Summary table creation")
summary_tab = annotation
for (cond in names(TAB)){
  tab = TAB[[cond]]
  mini_tab = cbind(tab$ID, tab$log2FoldChange, tab$padj, tab$REGULATION)
  colnames(mini_tab) = c("ID", paste(cond, c("log2FC", "padj", "REG"), sep = "_"))
  summary_tab = merge(summary_tab, mini_tab, by = "ID", all = T)
}
colnames(TurboPGM) = c("PROTEIN_NAME", paste0("TurboPGM_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGM, by = "PROTEIN_NAME", all = T)
colnames(TurboPGML4) = c("PROTEIN_NAME", paste0("TurboPGML4_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGML4, by = "PROTEIN_NAME", all = T)

write.table(summary_tab,paste0("Analyse/",file_name,"/",condition,"/Summary_",condition,".tab"), sep = "\t", row.names = F) 

### Venn Diagrame ####
print("Venn Digramm in progress")
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
LIST = c(UP_PKX, Intermediate_peak = list(AUTOGAMY$inter_peak))
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
            Intermediate_peak = AUTOGAMY$inter_peak)
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
LIST = c(UP_PKX, Early_peak = list(AUTOGAMY$early_peak))
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
            Early_peak = AUTOGAMY$early_peak)
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
LIST = list(DOWN_CTIP = stdCTIP[["DOWN_CTIP"]],
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
LIST = c(UP_PKX, DOWN_CTIP = list(stdCTIP[["DOWN_CTIP"]]))
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
            DOWN_CTIP = stdCTIP[["DOWN_CTIP"]],
            Intermediate_peak = AUTOGAMY$inter_peak)
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
            DOWN_CTIP = stdCTIP[["DOWN_CTIP"]],
            TurboPGM_PGML4 = turbo,
            Intermediate_peak = AUTOGAMY$inter_peak)
png(paste0(path,"Venn_CTIP_PKX_INTER_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avev les CTIP seuls
LIST = list(Early = downCTIP$CTIP_early,
            Intermediate = downCTIP$CTIP_inter)
png(paste0(path,"Venn_CTIP.png"))
ggvenn(LIST,
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7)
dev.off()

# Avec les CTIP UP dérégulés + early peak
LIST = list(UP_PKX = up_pkx,
            DOWN_CTIP = stdCTIP[["DOWN_CTIP"]],
            Early_peak = AUTOGAMY$early_peak)
png(paste0(path,"Venn_CTIP_PKX_EARLY.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Avec les CTIP UP dérégulés + early peak + Turbo
LIST = list(UP_PKX = up_pkx,
            DOWN_CTIP = stdCTIP[["DOWN_CTIP"]],
            TurboPGM_PGML4 = turbo,
            Early_peak = AUTOGAMY$early_peak)
png(paste0(path,"Venn_CTIP_PKX_EARLY_TURBO.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()


### Répartition des profils selon les filtres ####
print("Profiles repartition barplot")
path = paste0(save_path,"Barplot_profil/")
dir.create(path ,recursive=T,showWarnings=F)

# Sur UP PGM KU80c & XRCC4 
UP_PKX = c(UP_PKX, UP_ALL = list(up_pkx))

Profile_Barplot(UP_PKX, "UP", path)
Profile_EnrichmentBarplot(UP_PKX, path,"UP" )

# Sur DOWN CTIP + UP PKX 

Profile_Barplot(stdCTIP, "CTIP", path)
Profile_EnrichmentBarplot(stdCTIP, path, "CTIP")

###### Répartition des gènes avec IES selon les filtres ####
print("Gene with IES repartition barplot")
path = paste0(save_path,"Barplot_IES/")
dir.create(path ,recursive=T,showWarnings=F)

# Sur UP PGM KU80c & XRCC4
IES_Barplot(UP_PKX, path, "UP")
IES_EnrichmentBarplot(UP_PKX, path)


# Sur DOWN CTIP + UP PKX
IES_Barplot(stdCTIP, path, "CTIP")
IES_EnrichmentBarplot(stdCTIP, path)



### Analyse statitique ####
source("0_Stat_function.R")

LIST = c(UP_PKX, stdCTIP)
names(LIST)= c(names(UP_PKX), names(stdCTIP))

#### Enrichissement en gènes intermediate peak parmis les UP et CTIP

my_data = cbind(annotation$ID, annotation$EXPRESSION_PROFIL)

sink(paste0(save_path,"/Enrichissement_profile.txt"))
Enrichment_padj(LIST, my_data)
sink()


#### Imprimer l'état de R ####
sink(paste0(save_path,"/sessionInfo.txt"))
print(sessionInfo())
sink()