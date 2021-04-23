options(stringsAsFactors = FALSE)

annotation = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
annotation$SYNONYMS[grep("PTET.51.1.G0490126",annotation$ID)]="PGM"
annotation$SYNONYMS[grep("PTET.51.1.G0110267",annotation$ID)]="PGML1"
annotation$SYNONYMS[grep("PTET.51.1.G0380073",annotation$ID)]="PGML2"
annotation$SYNONYMS[grep("PTET.51.1.G0010374",annotation$ID)]="PGML3a"
annotation$SYNONYMS[grep("PTET.51.1.G0080308",annotation$ID)]="PGML3b"
annotation$SYNONYMS[grep("PTET.51.1.G0020217",annotation$ID)]="PGML3c"
annotation$SYNONYMS[grep("PTET.51.1.G0340197",annotation$ID)]="PGML4a"
annotation$SYNONYMS[grep("PTET.51.1.G0480099",annotation$ID)]="PGML4b"
annotation$SYNONYMS[grep("PTET.51.1.G0570051",annotation$ID)]="PGML5a"
annotation$SYNONYMS[grep("PTET.51.1.G0510172",annotation$ID)]="PGML5b"

seqlength=read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.transcript.fa.seqlength",h=T)
rownames(seqlength)=sub("PTET.51.1.T","PTET.51.1.G",seqlength$ID)

gene_autogamy=read.table("DATA/autogamy_ptetraurelia_mac_51_annotation_v2.0_significant.tab", header=T, sep="\t")
annotation = merge(annotation[,c(1,3)],gene_autogamy[,c(1,6)],by.x = "ID", by.y = "ID")

for (condition in list.files("./DATA/RNAseq") ){

  file_path=paste0("./DATA/RNAseq/", condition)
  files = list.files(path = file_path,pattern ="TOPHAT.pt_51.MIN_QUAL_30.htseq.count")
  
  tab = annotation
  rpm = annotation
  rpkm = annotation
  names = c("ID","SYNONYMS","EXPRESSION_PROFIL")
  for (i in files){
    a = read.table(paste0(file_path,"/",i))
    names = c(names,sub(".TOPHAT.pt_51.MIN_QUAL_30.htseq.count","",sub("PTET_mRNA_","",i )))
    tab = merge(tab,a, by.x = "ID", by.y = "V1")
    
    mapped_reads=sum(a$V2)
    a$RPM= a$V2 / mapped_reads *1e6
    a$RPKM= (a$V2 *1e3) / (seqlength[a$V1,]$LENGTH * (mapped_reads/1e6) )
    
    rpm = merge(rpm,a[,c(1,3)], by.x = "ID", by.y = "V1")
    rpkm = merge(rpkm,a[,c(1,4)], by.x = "ID", by.y = "V1")
    
  }
  colnames(tab)=names
  colnames(rpm)=names
  colnames(rpkm)=names
  
  
  write.table(tab,paste0("./count/",condition,"_expression_table.tab"), sep="\t",row.names=F,quote=F)
  write.table(rpm,paste0("./count_rpm/",condition,"_expression_table_RPM.tab"), sep="\t",row.names=F,quote=F)
  write.table(rpkm,paste0("./count_rpkm/",condition,"_expression_table_RPKM.tab"), sep="\t",row.names=F,quote=F)

}