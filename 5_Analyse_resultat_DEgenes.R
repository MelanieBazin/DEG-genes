####
# Analysis of the deregulated genes
# -> Expression profiles checking
# -> Venn diagram
# -> Proportion of the different profil cluster
####
options(stringsAsFactors = FALSE)
library(stringr) 
library(ggvenn)
library("VennDiagram")

source("0_Cluster.R") # Groups parameters
source("0_Functions.R") # Library of homemade function

# Access to the files folder
date = Sys.Date()
date = "2022-08-17"
condition =  names(rnai_list)[2]
path = paste0("./Analyse/",date,"_DESeq2_analysis/")
path = paste0(path,
              list.files(path)[grep(condition,list.files(path))],"/")


data_tab = read.table(paste0(path, condition,"_expression_table_vst.tab"))
data_mean = read.table(paste0(path, condition,"_MEANexpression_table_vst.tab"))

# Definition of the saving space
save_path = paste0(path, "Analyse/")
dir.create(save_path,recursive=T,showWarnings=F)

#### Creation of a summary table with all the data ####
print("Summary table creation")

# Merge table of deregulation and annotation information
summary_tab = annotation
for (cond in names(DEtab_list)){
  tab = DEtab_list[[cond]]
  mini_tab = cbind(tab$ID, tab$log2FoldChange, tab$padj, tab$REGULATION)
  colnames(mini_tab) = c("ID", paste(cond, c("log2FC", "padj", "REG"), sep = "_"))
  summary_tab = merge(summary_tab, mini_tab, by = "ID", all = T)
}
rm(mini_tab,cond)

# Add the information of turoID data (Guérinau et al.2020)
TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')

colnames(TurboPGM) = c("PROTEIN_NAME", paste0("TurboPGM_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGM, by = "PROTEIN_NAME", all = T)
colnames(TurboPGML4) = c("PROTEIN_NAME", paste0("TurboPGML4_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGML4, by = "PROTEIN_NAME", all = T)

write.table(summary_tab,paste0(path,"Summary_",condition,".tab"), sep = "\t", row.names = F)

#### Draw expression profile for a set of genes ####
selection = c("Ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", 
              "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")
selection = sort(selection)

select_ID =c()
name = c()
for( i in 1:length(selection)){
  select_ID = c(select_ID,annotation$ID[grep(selection[i],annotation$Name, ignore.case = T)]) 
  name = c(name,annotation$Name[grep(selection[i],annotation$Name, ignore.case = T)]) 
}
names(select_ID)=name
rm(selection,i,name)

ExpressionProfils(type = "vst",
                  condition = condition,
                  path = path,
                  name = "forRNAi",
                  select_ID = select_ID)

### Comparison of deregulated genes lists ####

source("5-1_Filters_candidats.R") # Comparison filter
UP_PKX.DOWN_C = Crossinglist(UP_filter,DOWN_filter)

##### Venn Diagram ####
print("Venn Diagramm in progress")
save_path2 = paste0(save_path,"VennDiagrame/")
dir.create(save_path2 ,recursive=T,showWarnings=F)

# Crossing LATE UP-deregulated genes in PGM, KU80c and XRCC4 time course
a1 = UP_filter$UP_KU80c
a2 = UP_filter$UP_PGM
a3 = UP_filter$UP_XRCC4

pdf(paste0(save_path2,"VennProp_UPpkx.pdf"))
grid.newpage() 
overrideTriple=T
draw.triple.venn(area1 = length(a1),
                 area2 = length(a2),
                 area3 = length(a3),
                 n12 = length(intersect(a1, a2)),
                 n23 = length(intersect(a2, a3)),
                 n13 = length(intersect(a1, a3)),
                 n123 = length(intersect(intersect(a1, a2),a3)),
                 category = c(paste("KU80c KD",length(a1), sep = "\n"),
                              paste("PGM KD",length(a2), sep = "\n"),
                              paste("XRCC4 KD",length(a3), sep = "\n")),
                 fill = c("deeppink1","dodgerblue3","ivory2"),
                 alpha = rep(0.5, 3),
                 scaled = TRUE)
dev.off()

# Crossing LATE UP-deregulated genes with DOWN-dergulated genes in CTIP time course
a1 = UP_filter$UP_PKX
a2 = DOWN_filter$DOWN_CTIP

pdf(paste0(save_path2,"VennProp_UPpkx_DOWNc.pdf"))
grid.newpage() 
draw.pairwise.venn(area1 = length(a1),
                   area2 = length(a2),
                   cross.area = length(intersect(a1, a2)),
                   category = c(paste("UPpkx",length(a1), sep = "\n"),
                                paste("DOWNc",length(a2), sep = "\n")),
                   fill = c("plum4", "chartreuse3"),
                   alpha = rep(0.5, 2),
                   scaled = TRUE)
dev.off()


##### Expression cluster representation depending on the filter ####
print("Profiles repartition barplot")
save_path2 = paste0(save_path,"Barplot_profil/")
dir.create(save_path2 ,recursive=T,showWarnings=F)

Profile_Barplot(UP_filter, "UP_PKX", save_path2, w= 14)
Profile_Barplot(DOWN_filter, "DOWN_C", save_path2)
Profile_Barplot(UP_PKX.DOWN_C, "UPpkx-DOWNc", save_path2, w= 40)
Profile_Barplot(DE_genes, "Candidates", save_path2)

# Statistic of the enrichment
my_data = cbind(annotation$ID, annotation$EXPRESSION_PROFIL)

sink(paste0(save_path2,"/Enrichissement_profile_UP_PKX.txt"))
Enrichment_padj(UP_filter, my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_DOWN_C.txt"))
Enrichment_padj(DOWN_filter, my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_Candidates.txt"))
Enrichment_padj(DE_genes, my_data)
sink()


##### IES+ genes among the filters ####
print("Gene with IES repartition barplot")
save_path2 = paste0(save_path,"Barplot_IES/")
dir.create(save_path2 ,recursive=T,showWarnings=F)

IES_Barplot(UP_filter, "UP_PKX", save_path2, w= 14)
IES_Barplot(DOWN_filter, "DOWN_C", save_path2)
IES_Barplot(UP_PKX.DOWN_C, "UPpkx-DOWNc", save_path2, w= 40)
IES_Barplot(DE_genes, "Candidates", save_path2)
IES_Barplot(AUTOGAMY, "ExpressionProfile", save_path2, w= 14)

# Statistic of the enrichment
my_data = cbind(annotation$ID, annotation$NB_IES)
my_data[my_data[,2] != 0,2] = "IES+"
my_data[my_data[,2] == 0,2] = "IES-"

sink(paste0(save_path2,"/Enrichissement_IES_UP_PKX.txt"))
Chi2_pvalue(UP_filter, my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_IES_DOWN_C.txt"))
Chi2_pvalue(DOWN_filter, my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_IES_Candidates.txt"))
Chi2_pvalue(DE_genes, my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_IES_ExpProfile.txt"))
Chi2_pvalue(DE_genes, my_data)
sink()

#### Print R status ####
sink(paste0(save_path,"/sessionInfo.txt"))
print(sessionInfo())
sink()