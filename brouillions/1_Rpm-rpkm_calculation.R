options(stringsAsFactors = FALSE)


seqlength=read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.transcript.fa.seqlength",h=T)
rownames(seqlength)=sub("PTET.51.1.T","PTET.51.1.G",seqlength$ID)


for (condition in list.files("./DATA/EXPRESSION") ){
  
  tab = read.table(paste0("./DATA/EXPRESSION/",condition), h=T)
  
  rpm = matrix(data = NA,nrow = nrow(tab),ncol = ncol(tab))
  colnames(rpm) = colnames(tab)
  rownames(rpm)=rownames(tab)
  
  rpkm = matrix(data = NA,nrow = nrow(tab),ncol = ncol(tab))
  colnames(rpkm) = colnames(tab)
  rownames(rpkm)=rownames(tab)
  
  for (i in 1:(ncol(tab))){
    mapped_reads=sum(tab[,i])
    rpm[,i]= tab[,i] / mapped_reads *1e6
    rpkm[,i]= (tab[,i] *1e3) / (seqlength[rownames(tab),]$LENGTH * (mapped_reads/1e6) )
    
  }  

  write.table(rpm,paste0("./DATA/RPM/",sub(".tab","",condition),"_expression_table_RPM.tab"), sep="\t",row.names=T,quote=F)
  write.table(rpkm,paste0("./DATA/RPKM/",sub(".tab","",condition),"_expression_table_RPKM.tab"), sep="\t",row.names=T,quote=F)
  
}
