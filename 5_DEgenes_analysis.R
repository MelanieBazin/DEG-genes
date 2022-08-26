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
condition =  names(rnai_list)[2]
path = paste0("./Analyse/",date,"_DESeq2_analysis/")
path = paste0(path,
              list.files(path)[grep(paste0(condition,"_FC-", FC, "_pval-", pvalue),list.files(path))],"/")


data_tab = read.table(paste0(path, condition,"_expression_table_vst.tab"))
data_mean = read.table(paste0(path, condition,"_MEANexpression_table_vst.tab"))

# Definition of the saving space
save_path = paste0(path, "Analyse/")
dir.create(save_path,recursive=T,showWarnings=F)

#### Draw expression profile for a set of genes ####
selection = c("NOWA","SPT", "PTIWI","DCL","mt","TFIIS4","Spo11", "CtIP" ,"Mre11",
              "CER","Rad51", "Ku","PGM", "XRCC4","Lig",
              "EZL","EED", "EAP","RF" ,"SUZ", "CAF1")
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
UP_PKX.DOWN_C = Crossinglist(UP_filter,list(DOWN_CTIP = DOWN_filter$DOWN_CTIP))
UP_PKX.DOWN_C = UP_PKX.DOWN_C[rapply(UP_PKX.DOWN_C, length, how="list") != 0]

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

# Crossing LATE UP-deregulated genes with DOWN-deregulated genes in CTIP time course
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

Profile_Barplot(list(UP = UP_filter$UP_KU80c, DOWN = DOWN_pkx_filter$DOWN_KU80c), "KU80c", save_path2)
Profile_Barplot(list(UP = UP_filter$UP_PGM, DOWN = DOWN_pkx_filter$DOWN_PGM), "PGM", save_path2)
Profile_Barplot(list(UP = UP_filter$UP_XRCC4, DOWN = DOWN_pkx_filter$DOWN_XRCC4), "XRCC4", save_path2)
Profile_Barplot(list(UP = UP_C_filter$UP_CTIP, DOWN = DOWN_filter$DOWN_CTIP), "CTIP", save_path2)
Profile_Barplot(UP_PKX.DOWN_C, "UPpkx-DOWNc", save_path2, w= 30)
Profile_Barplot(DE_genes, "Candidates", save_path2)

# Statistic of the enrichment
my_data = cbind(annotation$ID, annotation$EXPRESSION_PROFIL)

sink(paste0(save_path2,"/Enrichissement_profile_KU80c.txt"))
Enrichment_padj(list(UP = UP_filter$UP_KU80c, DOWN = DOWN_pkx_filter$DOWN_KU80c), my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_PGM.txt"))
Enrichment_padj(list(UP = UP_filter$UP_PGM, DOWN = DOWN_pkx_filter$DOWN_PGM), my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_XRCC4.txt"))
Enrichment_padj(list(UP = UP_filter$UP_XRCC4, DOWN = DOWN_pkx_filter$DOWN_XRCC4), my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_CTIP.txt"))
Enrichment_padj(list(UP = UP_C_filter$UP_CTIP, DOWN = DOWN_filter$DOWN_CTIP), my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_UPpkx_DOWNc.txt"))
Enrichment_padj(UP_PKX.DOWN_C, my_data)
sink()
sink(paste0(save_path2,"/Enrichissement_profile_Candidates.txt"))
Enrichment_padj(DE_genes, my_data)
sink()


##### IES+ genes among the filters ####
print("IES repartition barplot")
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

sink(paste0(save_path2,"/Chi2_IES_UP_PKX.txt"))
Chi2_pvalue(UP_filter, my_data)
sink()
sink(paste0(save_path2,"/Chi2_IES_DOWN_C.txt"))
Chi2_pvalue(DOWN_filter, my_data)
sink()
sink(paste0(save_path2,"/Chi2_IES_Candidates.txt"))
Chi2_pvalue(DE_genes, my_data)
sink()
sink(paste0(save_path2,"/Chi2_IES_ExpProfile.txt"))
Chi2_pvalue(AUTOGAMY, my_data)
sink()

#### GO therm analysis ####
print("GO analysis")
save_path2 = paste0(path,"GO_therm/")
dir.create(save_path2 ,recursive=T,showWarnings=F)

interestGenes = candidats
analysis_name = "Candidats"

##### Plot the different therms ####
gene_tab = "./DATA/GOenrichment/paramecium_gene_association.fb"
gene_tab=read.table(gene_tab,sep="\t",quote='',comment.char='')[,c(3,2,9,5,10)]
colnames(gene_tab) = c("ID","NAME","VOCABULARY","GO_ID","DESC")

gene_tab_slim = "./DATA/GOenrichment/paramecium_gene_association.slim.fb"
gene_tab_slim=read.table(gene_tab_slim,sep="\t",quote='',comment.char='')[,c(3,2,9,5,10)]
colnames(gene_tab_slim) = c("ID","NAME","VOCABULARY","GO_ID","DESC")

for(voc in c("MF","BP")){
  # Selection of the GO therm type
  if (voc == "MF"){
    tab = gene_tab[gene_tab$VOCABULARY == "F",]
    tab_slim = gene_tab_slim[gene_tab_slim$VOCABULARY == "F",]
  }else if (voc == "BP"){
    tab = gene_tab[gene_tab$VOCABULARY == "P",]
    tab_slim = gene_tab_slim[gene_tab_slim$VOCABULARY == "P",]
  }
  
  # Selection of the gene of interest
  tab.select = tab[which(is.element(tab$ID, interestGenes)),]
  tab.select = as.data.frame(table(tab.select$DESC))
  colnames(tab.select) = c("DESC","Nb")
  
  tab = as.data.frame(table(tab$DESC))
  colnames(tab) = c("DESC","Nb")

  # Table with the nb of genes in each category for all or selected genes
  GO_tab = merge(tab,tab.select, by = "DESC", all = T)
  colnames(GO_tab) = c("DESC", "ALL", analysis_name)
  write.table(GO_tab, paste0(save_path2,voc,"_",analysis_name,"_Nb_GOtherm.tab"),row.names=F,quote=F,sep="\t")
  
  # Barplot of the therm frequency for the selected genes
  tab.select$class = with(tab.select, reorder(DESC, Nb))
  gp = ggplot(tab.select, aes(class,Nb)) + 
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_classic()
  ggsave(paste0(save_path2,voc,"_",analysis_name,"_GOtherm.pdf"),
         height = 7, width = 8, device = "pdf", )
  
  ## Same with slim therms
  # Selection of the gene of interest
  tab_slim.select = tab_slim[which(is.element(tab_slim$ID, interestGenes)),]
  tab_slim.select = as.data.frame(table(tab_slim.select$DESC))
  colnames(tab_slim.select) = c("DESC","Nb")
  
  tab_slim = as.data.frame(table(tab_slim$DESC))
  colnames(tab_slim) = c("DESC","Nb")
  
  # Table with the nb of genes in each category for all or selected genes
  GO_tab_slim = merge(tab_slim,tab_slim.select, by = "DESC", all = T)
  colnames(GO_tab_slim) = c("DESC", "ALL", analysis_name)
  write.table(GO_tab, paste0(save_path2,voc,"_",analysis_name,"_Nb_GOtherm_slim.tab"),row.names=F,quote=F,sep="\t")
  
  # Barplot of the therm frequency for the selected genes
  tab_slim.select$class = with(tab_slim.select, reorder(DESC, Nb))
  gp = ggplot(tab_slim.select, aes(class,Nb)) + 
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_classic()
  ggsave(paste0(save_path2,voc,"_",analysis_name,"_GOtherm_slim.pdf"),
         height = 7, width = 8, device = "pdf", )
}

##### Identified enriched therms ####
print("GO analysis")
map_file_path = "./DATA/GOenrichment/ptetraurelia_mac_51_annotation_v2.0.protein.InterProScan.tab.gene_association.map"
map_file = readMappings(file =map_file_path)
names(map_file) = sub(".P",".G", names(map_file))

map_slim_path = "./DATA/GOenrichment/paramecium_gene_association.slim.fb"
gene_tab=read.table(map_slim_path,sep="\t",quote='',comment.char='')[,c(3,9,5)]
colnames(gene_tab) = c("ID","VOCABULARY","GO_ID")

for(voc in c("MF","BP")) {
  # Test with all GO therm
  enrichedGenes = top_go_enrichment(map_file,
                                    interestGenes,
                                    voc,
                                    pvalue=0.05, 
                                    adjust_pvalue=T,
                                    prefix=paste0(save_path2,voc,"_",analysis_name,"_"))
  write.table(enrichedGenes, paste0(save_path2,voc,"_",analysis_name,"_enriched_GO.tab"),row.names=F,quote=F,sep="\t")
  
  # Test with slim GO therm
  if (voc == "MF"){
    tab = gene_tab[gene_tab$VOCABULARY == "F",]
  }else if (voc == "BP"){
    tab = gene_tab[gene_tab$VOCABULARY == "P",]
  }
  map_file_slim = as.list(tab$GO_ID)
  names(map_file_slim) = tab$ID
  
  enrichedGenes = top_go_enrichment(map_file_slim,
                                    interestGenes,
                                    voc,
                                    pvalue=0.05, 
                                    adjust_pvalue=T,
                                    prefix=paste0(save_path2,voc,"_",analysis_name,"_",))
  write.table(enrichedGenes, paste0(save_path2,voc,"_",analysis_name,"_enriched_GOslim.tab"),row.names=F,quote=F,sep="\t")
}


#### Creation of a summary table with all the data ####
print("Summary table")
# Merge table of deregulation and annotation information
summary_tab = annotation
for (cond in names(DEtab_list)){
  tab = DEtab_list[[cond]]
  mini_tab = cbind(tab$ID, tab$log2FoldChange, tab$padj, tab$REGULATION)
  colnames(mini_tab) = c("ID", paste(cond, c("log2FC", "padj", "REG"), sep = "_"))
  summary_tab = merge(summary_tab, mini_tab, by = "ID", all = T)
}
rm(mini_tab,cond)

write.table(summary_tab,paste0(path,"Summary_",condition,".tab"), sep = "\t", row.names = F)

# Extracted the raw corresponding to the candidates
candidats_tab = summary_tab[which(is.element(summary_tab$ID, candidats)),]
write.table(summary_tab,paste0(path,"Summary-candidats_",condition,".tab"), sep = "\t", row.names = F)

# Add the information of turoID data (Guérinau et al.2020)
TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')

colnames(TurboPGM) = c("PROTEIN_NAME", paste0("TurboPGM_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGM, by.x = "Protein.name",  by.y = "PROTEIN_NAME", all = T)
colnames(TurboPGML4) = c("PROTEIN_NAME", paste0("TurboPGML4_", c("log2FC", "-log10pval")))
summary_tab = merge(summary_tab, TurboPGML4, by.x = "Protein.name", by.y = "PROTEIN_NAME", all = T)

write.table(summary_tab,paste0(path,"Summary-turbo_",condition,".tab"), sep = "\t", row.names = F)

#### Print R status ####
sink(paste0(save_path,"/sessionInfo.txt"))
print(sessionInfo())
sink()

