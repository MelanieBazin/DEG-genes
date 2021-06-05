library(stringr)
library(gplots)
source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")

data_path = "./DATA/RNAseq/"
save_path = "./Analyse/RowProfils/"

dir.create("./DATA/EXPRESSION/",recursive=T,showWarnings=F)
dir.create(save_path,recursive=T,showWarnings=F)

data_directories = list.files(data_path)

for (i in data_directories){
  list= list.files(paste0(data_path,i))
  
  # Extraire les timing des nom des fichier
  timing = str_split(list, '_', simplify=TRUE)
  for(a in c(1,4)){
   if(length(unique(timing[,a]))>1){
     time = a
   }
  }
  timing = timing[,time]
  x = str_split(timing, '-', simplify=TRUE)[1,1]
  timing = gsub(paste0(x,"-"), "",timing)
  
  # Ouvrir les fichier et faire un tableau par condition
  table = data.frame(annotation$ID)
  colnames(table)="ID"
  for(j in list){
    tab = read.table(paste0(data_path,i,"/",j),)
    colnames(tab)[1]="ID"
    table = merge(table, tab, by = "ID")
  }
  rownames(table)=as.character(table$ID)
  if (colnames(table)[1]=="ID"){
    table = table[,-1]
  }
  colnames(table) = timing
  # Réordonner les colonnes par ordre chonologique
  colorder = c(ncol(table),order(as.numeric(gsub("T", "", colnames(table)[-ncol(table)]))))
  table = table[,colorder]
  colnames(table)[1] = "Veg"
  
  # Retirer les lignes qui ne correspondent pas à un gènes
  table = table[which(is.element(rownames(table),annotation$ID)),]
  
  # Créer un nouveau tableau
  write.table(table,paste0("./DATA/EXPRESSION/",i,".tab"),sep="\t", row.names=T,quote=F)
} 
  # Faire plot, box, plot et heatmap sur les row data pour verifier s'il ya un problème
  if(!is.element(F, rownames(table)==annotation$ID)){
    data_log =  as.matrix(log(table+1))
    expr_profil = unique(annotation$EXPRESSION_PROFIL)
    
    for(p in expr_profil){
      id = annotation$ID[grep(p, annotation$EXPRESSION_PROFIL)]
      png(paste0(save_path,i,"_",p,"_Profil.png"))
        graph = plotGenes(table[id,], title = p, yMax = max(table[id,]))
        print(graph)
      dev.off()
      
      png(paste0(save_path,i,"_",p,"_Boxplot.png"))
        graph = boxplot(data_log[id,],main = p ,ylab = "log(EXPRESSION)",
                        xlab = colnames(data_log))
        print(graph)
      dev.off()
      
      png(paste0(save_path,i,"_",p,"_Heatmap.png"))
        graph = heatmap.2(data_log[id,], Colv = NULL, trace = 'none', dendrogram = 'none')
        print(graph)
      dev.off()
    }

  }
}


