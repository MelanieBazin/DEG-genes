options(stringsAsFactors = FALSE)

OpenDataCount <- function(path, condition){
  files = list.files(path = path)
  
  if (condition == "PGM" | condition == "KU80c") { 
    control_nd7 = read.table(paste0(path,"/",files[grep("ND7",files)]),h = T,sep = "\t")
    control_nd7 = control_nd7[,c(1,ncol(control_nd7),2:(ncol(control_nd7)-1))]
    colnames(control_nd7) = paste(colnames(control_nd7),"CTRL_ND7", sep = "_" )
    
    control_icl7 = read.table(paste0(path,"/",files[grep("ICL7",files)]),h = T,sep = "\t")
    control_icl7 = control_icl7[,c(1,ncol(control_icl7),2:(ncol(control_icl7)-1))]
    colnames(control_icl7) = paste(colnames(control_icl7),"CTRL_ICL7", sep = "_" )
    
    control = merge(control_nd7, control_icl7, by.x = "ID_CTRL_ND7", by.y = "ID_CTRL_ICL7")
    
    rnai = read.table(paste0(path,"/",files[grep(condition,files)]),h = T,sep = "\t")
    rnai = rnai[,c(1,ncol(rnai),2:(ncol(rnai)-1))]
    colnames(rnai) = paste(colnames(rnai),"RNAi", sep = "_" )
    
    countdata=merge(control, rnai, by.x = "ID_CTRL_ND7", by.y = "ID_RNAi")
    row.names(countdata) = countdata$ID_CTRL_ND7
    countdata = countdata[,-1]
    
  } else if (condition == "XRCC4" | condition == "CTIP"){
    files = files[grep(condition,files)]
    control = read.table(paste0(path,"/",files[grep("CTRL",files)]),h = T,sep = "\t")
    control = control[,c(1,ncol(control),2:(ncol(control)-1))]
    colnames(control) = paste(colnames(control),"CTRL", sep = "_" )
    
    rnai = read.table(paste0(path,"/",files[setdiff(1:length(files),grep("CTRL",files))]),h = T,sep = "\t")
    rnai = rnai[,c(1,ncol(rnai),2:(ncol(rnai)-1))]
    colnames(rnai) = paste(colnames(rnai),"RNAi", sep = "_" )
    
    countdata=merge(control, rnai, by.x = "ID_CTRL", by.y = "ID_RNAi")
    row.names(countdata) = countdata$ID_CTRL
    countdata = countdata[,-1]
    
  
   } else {
    print(paste("Condition",condition,"non pris en charge"))
   }
  
  return(countdata)
}


CreationInfoData <- function(countdata){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 3)
  
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Noms", "Feeding", "Timing")
  
  
  timing = sub("RNAi","",colnames(countdata))
  timing = sub("CTRL","",timing)
  timing = sub("ND7","",timing)
  timing = sub("ICL7","",timing)
  timing = gsub("_","",timing)
  
  if (condition == "PGM"|condition=="KU80c"){
    infodata[,"Feeding"] = c(rep("ND7 PGM/KU80c",length(grep("ND7", colnames(countdata)))),
                             rep("ICL7 PGM/KU80c",length(grep("ICL7", colnames(countdata)))),
                             rep(paste("RNAi",condition),length(grep("RNAi", colnames(countdata))))) 
  }else{
    infodata[,"Feeding"] = c(rep(paste("ND7", condition),length(grep("CTRL", colnames(countdata)))),
                             rep(paste("RNAi",condition),length(grep("RNAi", colnames(countdata)))))
  }
  
  
  infodata[,"Noms"] = colnames(countdata)
  infodata[,"Timing"] = timing
  
  return(infodata)
}

