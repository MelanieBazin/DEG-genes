library("RColorBrewer")
library("vsn")
library("gplots")
library("pheatmap")
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)


library("vsn")

library(seqinr)
hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(255)
#hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol = rev(hmcol)

options(stringsAsFactors = FALSE)

library(gplots)
library(pheatmap)

f<-function(x) { x }

library(latticeExtra)

library(seqinr)
#############################################
# PCA
#############################################


library(FactoMineR)
library(gtools)

require(ggplot2)

PCA_plot_generator <- function(Expression_Mat, colors,max_dim=3,barplot_max_dim=3,image_prefix="PCA_",show_barplot=T, vline=0, ...) {
  resExp = PCA(t(Expression_Mat), graph = F)
  if(show_barplot) {
    eigenvalues <- resExp$eig
    pdf(paste0(image_prefix,"_PCA_Variance.pdf"))
    barplot(eigenvalues[1:barplot_max_dim, 2], names.arg=1:barplot_max_dim, 
            main = "Variances",
            xlab = "Principal Components",
            ylab = "Percentage of variances",
            col ="gray", ...)
    if(vline!=0) {
      abline(v=vline,lty=2,lwd=2)
    }
    dev.off()
  }    
  
  for (i in 1:dim(combn(1:max_dim,2))[2]) {
    
    gp<-plot.PCA(resExp, axes = combn(1:max_dim,2)[,i], habillage = "ind", col.hab = colors, ...)
    ggsave(paste0(image_prefix,i,".pdf"), plot = gp)
  }
  
  
}


# cut a vector into bins 
evenbins <- function(x, bin.count=10, order=T) {
  bin.size <- rep(length(x) %/% bin.count, bin.count)
  bin.size <- bin.size + ifelse(1:bin.count <= length(x) %% bin.count, 1, 0)
  bin <- rep(1:bin.count, bin.size)
  if(order) {    
    bin <- bin[rank(x,ties.method="random")]
  }
  return(factor(bin, levels=1:bin.count, ordered=order))
}


# adhoc function to get autogmay gene enrichment
get_autogamy_enrichment<-function(selected_ids) {
  
  wt_autogamy_selected=wt_autogamy[selected_ids,]
  induced_profiles = profiles[grep("repression",profiles,invert=T)]
  induced_profiles = induced_profiles[induced_profiles!="none"]
  line=c(
    # Number
    length(selected_ids),
    # Autogamy genes
    dim(wt_autogamy_selected[wt_autogamy_selected$SIGNIFICANT,])[1],
    round(dim(wt_autogamy_selected[wt_autogamy_selected$SIGNIFICANT,])[1]/length(selected_ids)*100,2),
    
    # induced
    sum(table(wt_autogamy_selected$EXPRESSION_PROFIL)[induced_profiles]),
    round(sum(table(wt_autogamy_selected$EXPRESSION_PROFIL)[induced_profiles])/dim(wt_autogamy_selected[wt_autogamy_selected$SIGNIFICANT,])[1]*100,2),      
    # repressed
    sum(table(wt_autogamy_selected$EXPRESSION_PROFIL)[profiles[grep("repression",profiles)]]),
    round(sum(table(wt_autogamy_selected$EXPRESSION_PROFIL)[profiles[grep("repression",profiles)]])/dim(wt_autogamy_selected[wt_autogamy_selected$SIGNIFICANT,])[1]*100,2),      
    # by groups
    as.vector(table(wt_autogamy_selected$EXPRESSION_PROFIL)[profiles] ),
    round(as.vector(table(wt_autogamy_selected$EXPRESSION_PROFIL)[profiles] )/dim(wt_autogamy_selected[wt_autogamy_selected$SIGNIFICANT,])[1]*100,2)
    
  )
  line
}

# internal chi2 function
chi2 <- function(v1, v2, correct=F){
  
  m = matrix(data=NA, ncol=2, nrow=2)
  #print(paste(v1[1],"/",v1[2], v2[1],"/",v2[2],sep=" "))
  m[ , 1] = c(v1[1], v2[1])
  m[ , 2] = c( v1[2]-v1[1] , v2[2]-v2[1] )
  
  chi2 = chisq.test(m, correct=correct)
  chi2
}

# my chi square test
my_chi2<-function(x,ctl) {
  t=chi2(as.numeric(c(ctl[1],ctl[2])),as.numeric( c(x[1],x[2])))
  pv=t$p.value
  signif=""
  if(pv < 1e-200) {
    signif="****"
  } else {
    if(pv < 1e-100) {
      signif="***"
    } else {   
      
      if(pv < 1e-20) {
        signif="**"
      } else {
        if(pv < 1e-10) {
          signif="*"
        }
      }
    }
  }
  print(paste(x[1],x[2],paste0(round(x[1]/x[2]*100,1),"%"),ctl[1],ctl[2],format(pv,digits=3),format(t$statistic,digits=1),signif,sep=" "))
}




library(topGO)

# GO enrichment analysis using library topGO
top_go_enrichment <- function(map_file,myInterestingGenes,vocabulary="MF", random=FALSE, pvalue=0.05, adjust_pvalue=TRUE,prefix="", geneUniverse=c()) {
  
  
  
  geneID2GO <- readMappings(file =map_file)
  
  #str(head(geneID2GO))
  GO2geneID <- inverseList(geneID2GO)
  #str(head(GO2geneID))
  
  geneNames <- names(geneID2GO)
  if(random==TRUE) {
    print("Random")
    myInterestingGenes <- sample(geneNames, length(myInterestingGenes))
  }
  
  if(length(geneUniverse) != 0) {
    myInterestingGenes = intersect(geneNames,geneUniverse)
    
  }
  
  #head(geneNames)
  
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  #str(geneList)
  GOdata <- new("topGOdata", ontology = vocabulary, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  
  results.table <- GenTable(GOdata, classicFisher = resultFisher , topNodes = length(resultFisher@score))
  
  results.table.bh = results.table[results.table$classicFisher<=pvalue,]
  if(adjust_pvalue==TRUE) {
    results.table.bh = results.table[which(p.adjust(results.table[,"classicFisher"],method="BH")<=pvalue),] 
  }
  
  
  results.table.bh
  #~    write.table(results.table.bh,paste(prefix,"top_go_enrichment_",basename(map_file),"_",vocabulary,".tab",sep=""),row.names=F,quote=F,sep="\t")
  
  
  
  #~    #We can also look at multiple GO terms at the same time:
  #~    if(dim(results.table.bh)[1]!=0) {
  #~    GOids.of.interest = results.table.bh[,"GO.ID"]
  #~    all.term.genes = genesInTerm(GOdata, GOids.of.interest)
  #~    # Which of these genes is in the bicluster?
  #~    genes.of.interest <- sapply(names(all.term.genes),function(x){intersect(all.term.genes[[x]],myInterestingGenes)})
  #~    #write.table(unlist(lapply(genes.of.interest, paste, collapse=" ")),"t",quote=F,sep="\t",col.names=F)
  #~    unique(as.vector(unlist(genes.of.interest)))
  #~    }
}


library(ggrepel)
library(viridis)
library(ggplot2)
# usefull function to get density points, need this function for ggplots2 plots
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
