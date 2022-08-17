options(stringsAsFactors = FALSE)
library(ggvenn)
library(ggplot2) 
library(RColorBrewer)

path = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/analyseDE/NewGraphs_UP_CTIP_inter/"
path = "./NewGraphs/"
dir.create(path ,recursive=T,showWarnings=F)


file = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/analyseDE/Motif/From_TSS_IN_MAC/UP_CTIP_inter/FIMO_1E-4/p-value_1.4e-05/Summary2_analyseDE.tab"
file = "./Summary2_analyseDE.tab"

summary = read.table(file, header = T, sep = "\t")

KU = summary$ID[which(summary$KU80c_REG == "Up-regulated")]
PGM = summary$ID[which(summary$PGM_REG == "Up-regulated")]
XRCC4 = summary$ID[which(summary$XRCC4_REG == "Up-regulated")]
CTIP_early = summary$ID[which(summary$CTIP_early_REG == "Down-regulated")]
CTIP_inter = summary$ID[which(summary$CTIP_inter_REG == "Down-regulated")]

UP = intersect(intersect(KU, PGM), XRCC4)
CTIP = unique(c(CTIP_early, CTIP_inter))

UP_CTIP = intersect(UP, CTIP)

INTER = summary$ID[which(summary$EXPRESSION_PROFIL == "Intermediate peak")]

UP_CTIP_inter = intersect(UP_CTIP, INTER)
write(UP_CTIP_inter, paste0(path,"ID_UPpkx_DOWNc_inter.txt"))

source("0_Visualisation_fonction.R")

LIST = list(
  PGM  = PGM,
  KU80c = KU,
  XRCC4 = XRCC4,
  CTIP = CTIP,
  UP_CTIP = UP_CTIP)

Profile_Barplot(LIST, "", path)

source("0_Stat_function.R")
my_data = as.data.frame(cbind(summary$ID, summary$EXPRESSION_PROFIL))

sink(paste0(path,"/Enrichissement_position_autogamy.txt"))
Enrichment_padj(LIST, my_data)
sink()
