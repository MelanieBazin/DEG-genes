###################
#### DESeq2 result extraction
#### -> Extract the comparison of expression level of a KD compared to a control for each autogamy stages
#### -> Visualization of the data by Volcanoplot
###################
options(stringsAsFactors = FALSE)

# For each autogamy stage
for(i in names(stages)) {
  print(paste(RNAi, "-",i, "stage"))
  
  c1=stages[[i]][grep("ctrl",stages[[i]])]
  c2=stages[[i]][grep(RNAi,stages[[i]])]
  if(length(c2)==0|length(c1)==0){}else{
    
    # Comparison of a given time course to the controls at a given autogamy stage
    resContrast=results(deseq,contrast=c("Condition",c2, c1))
    resContrast=resContrast[rowSums(counts(deseq)) > 0,] # Suppress the genes with no coverage
    
    # Define genes that reach the significant FC and p-value (define earlier)
    res = resContrast[!is.na(resContrast$padj) & resContrast$padj < pvalue ,]
    res = res[res$log2FoldChange >= log2(FC) | 
                res$log2FoldChange <= log2(1/FC),]
    
    res = as.data.frame(res)
    res=res[,c("baseMean","log2FoldChange","padj")]
    res$REGULATION=ifelse(res$log2FoldChange>0,"Up-regulated","Down-regulated")
    
    # Create an annotated version of the results table
    res.annot= merge(res,annotation,by.x="row.names",by.y="ID")
    res.annot=res.annot[order(res.annot$padj), ]
    colnames(res.annot)[1]="ID"
    write.table(res.annot,paste0(res_dir,"DEgenes_",condition,"_",i,".tab"),sep="\t",quote=F,row.names=F)

    #### Generation of Volcanoplots ####
    if(dim(res)[1] !=0) {
      res_vp=as.data.frame(resContrast)
      res_vp$SIGNIFICANT=FALSE
      res_vp[which(is.element(rownames(res_vp),rownames(res))),]$SIGNIFICANT=TRUE
      
      # Volcanoplot
      pdf(paste0(img_dir,"Volcanoplot_",condition,"_",i,".pdf"), width = 6, height = 6)
      plot(res_vp$log2FoldChange,-log(res_vp$padj),
           log="y",
           col=ifelse(res_vp$SIGNIFICANT,"indianred","gray"),
           xlab=paste0("log2(",c1,"/",c2,")"),
           ylab="-log(p-value)",
           pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      dev.off()
      
      # Table corresponding to Volcanoplot data
      all = length(annotation$ID)
      volca_tab = table(res$REGULATION)
      volca_tab = rbind(volca_tab,volca_tab/all*100)
      write.table(volca_tab, paste0(img_dir,"Volcanoplot_tab_",condition,"_",i,".tab"),sep="\t",quote=F,row.names=F)
      
      
    }
  }  
}
