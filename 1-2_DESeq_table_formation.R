####
# Put all data that need will be analysed together in a table
# + corrected the expression count by ComBat-seq
# + creation of the infodata table containing all information about the samples
####

options(stringsAsFactors = FALSE)
library(sva)
library(stringr)
library(gplots)

# set.seed(10111)

source("0_Functions.R") # Library of homemade function
annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")
profile_color = c("Early peak" = "purple3",
                  "Intermediate peak" = "red2",
                  "Late peak" = "chartreuse4",
                  "Early repression" = "dodgerblue3",
                  "Late induction" = "deeppink",
                  "Late repression" = "darkorange",
                  "none" = "snow3")

#### Creation of expression files form the count ####
# Definition of the identifier for each time points
RNAi_count_files = list(ICL7 = "ICL7-",
                        ICL7bis = c("DUHA140","DUHA141","DUHA142","DUHA143"),
                        EZL1 = "EZL1-",
                        EZL1bis = c("DUHA144","DUHA145","DUHA146","DUHA147"),
                        ND7_L = c("BET32","BET33","BET34", "BET35", "BET87","BET36", "BET37"),
                        ND7_K = "ND7-",
                        ND7_X = c("BET38","BET39","BET40", "BET41", "BET42","BET86"),
                        PGM = "PGM-",
                        KU80c = "KU80-",
                        XRCC4 = c("BET27","BET28","BET29", "BET30", "BET31","BET85"),
                        CTIP = c("BET22","BET23", "BET24","BET25","BET26"),
                        LIG4 = "BET88")
RNAi_count_files = RNAi_count_files[-grep("LIG",names(RNAi_count_files))]


data_path = "/Plateformes/PARAMECIUM/DATA/PARAMECIUM/COVERAGE/tetraurelia/ptetraurelia_mac_51/ptetraurelia_mac_51_annotation_v2.0/mRNA/htseq/"
data_files= list.files(data_path)
data_files = data_files[-grep("REMAP", data_files)]
data_files = data_files[-grep("count2", data_files)]
data_files = data_files[grep("TOPHAT", data_files)]

save_path = "./DATA/EXPRESSION/Profils/"
dir.create(save_path,recursive=T,showWarnings=F)

for (i in names(RNAi_count_files)){
  print(paste("Expression table for",i))
  # Extract the file names
  list = c()
  bar_code = RNAi_count_files[[i]]
  for (bc in bar_code){
    list = c(list,
             data_files[grep(bc,data_files)])
  }
  
  # Extract the timing form file names
  timing = str_split(list, '_', simplify=TRUE)
  for(a in c(4,2,1)){
    if(length(unique(timing[,a]))>2){
      time = a
    }
  }
  
  # Open files and create a table per KD
  table = data.frame(annotation$ID)
  colnames(table)="ID"
  for(j in list){
    tab = read.table(paste0(data_path,j),sep = '\t')
    tab = tab[grep("PTET", tab$V1),] # Remove lines not corresponding to gene
    
    t = str_split(j, '_', simplify=TRUE)[time]
    x = str_split(t, '-', simplify=TRUE)[1,1]
    t = gsub(paste0(x,"-"), "",t)
    
    colnames(tab)=c("ID",t)
    table = merge(table, tab, by = "ID")
  }
  # Put gene ID in row names
  rownames(table)=as.character(table$ID)
  if (colnames(table)[1]=="ID"){
    table = table[,-1]
  }
  
  # Reorder column by chronology (Veg in first)
  veg = grep("v",colnames(table),ignore.case = T)
  if (length(veg) > 0){if(veg == ncol(table)){
    colorder = c(ncol(table),order(as.numeric(gsub("T", "", colnames(table)[-ncol(table)]))))
    table = table[,colorder]
  }}else {
    colorder = c(order(as.numeric(gsub("T", "", colnames(table)))))
    table = table[,colorder]
  }
  
  colnames(table)[grep("v",colnames(table),ignore.case = T)] = "Veg"
  
  # Remove lines not corresponding to annotated genes
  table = table[which(is.element(rownames(table),annotation$ID)),]
  
  print(colnames(table))
  
  # Save the table
  write.table(table,paste0("./DATA/EXPRESSION/",i,".tab"),sep="\t", row.names=T,quote=F)
  
  # Plot all data on boxplot
  profils = unique(annotation$EXPRESSION_PROFIL)
  
  pdf(paste0(save_path, i, "_Boxplot_Expression.pdf"), width = 20, height = 8)
  par(mfrow = c(2,round(length(profils)/2)))
  for (p in profils){
    
    data_box = table[annotation$ID[annotation$EXPRESSION_PROFIL == p],]
    b = boxplot(data_box, plot = F)$stats
    mediane = b[3,]
    
    plot(NULL, main = p, sub = i,
         ylim = c(0, round(max(b[5,]), digits = -0.01)),
         ylab = "Expression level",
         xlim = c(0.5,ncol(b)+0.5),
         xlab = "",
         xaxt = "n")
    boxplot(data_box, add = T, boxwex =0.5, col = "white", outline = F,
            names = str_remove_all(colnames(data_box), paste0("_",i)))
    lines(mediane, col = profile_color[p], lwd = 4)
  }
  dev.off()
  
}

#### Setting files for DESeq2 analysis + batch correction ####
source("0_Cluster.R") # Groups parameters

# Creation of saving folder
dir.create("./DATA/For_DESeq/",recursive=T,showWarnings=F)

for (i in names(rnai_list)){
  # Opening files
  countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[i]])
  
  # Before ComBat-seq correction
  countdata = as.matrix(countdata)
  write.table(countdata,paste0("./DATA/For_DESeq/",i,"_expression_table_uncorrected.tab"), sep="\t",row.names=T,quote=F)
  
  # After ComBat-seq correction
  print(paste("Batch correction for",i))
  infodata = CreatInfoData(countdata, conditions = i , rnai_list, cluster)
  batch = paste(infodata$Seq_method,infodata$Labo, sep = "_")
  countdata = ComBat_seq(countdata, batch = batch)
  
  write.table(countdata,paste0("./DATA/For_DESeq/",i,"_expression_table_corrected.tab"), sep="\t",row.names=T,quote=F)
}

#### Print R status ####
sink(paste0("./DATA/For_DESeq/sessionInfo.txt"))
print(sessionInfo())
sink()
