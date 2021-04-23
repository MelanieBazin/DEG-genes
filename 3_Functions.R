options(stringsAsFactors = FALSE)

##### Chragement des libraies necessaire à l'analyse ####
library(FactoMineR)
library(factoextra)
library(DESeq2)
library(gplots)
library(seqinr)
library("ggvenn")
library("RColorBrewer")
library("pheatmap")
#####

# Variable pour l'enrichissement en gène de l'autogamie
wt_autogamy =read.table("./DATA/autogamy_ptetraurelia_mac_51_annotation_v2.0_significant.tab",h=T,sep="\t")
rownames(wt_autogamy)=wt_autogamy$ID

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


