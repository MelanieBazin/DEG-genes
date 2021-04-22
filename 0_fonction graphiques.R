options(stringsAsFactors = FALSE)

annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

##### Faire les courbes #####

BoxDrawer = function(tab, type, files, parameter, zoom_max, pas){
  condition = sub(paste0("_expression_table_",type,".tab"),"",files)
  tab = na.omit(tab)
  temp = merge(annotation[,c(1,3)],tab, by.x = "ID", by.y = "ID")

  early = tab[temp$EXPRESSION_PROFIL=="Early peak",]
  intermediate = tab[temp$EXPRESSION_PROFIL=="Intermediate peak",]
  
  
  png(paste0("./Graph/",condition,"/Box_",condition,"_",type,"_",parameter,".png"),width = 2000, height = 2000)
  par(mfrow=c(2,2)) 
  boxplot(early[,-1],
          ylab="Expression UA",
          xlab="Temps autogamie",
          main=paste0("Expression early peak in ", condition),
          ylim=c(0,ceiling(max(apply(early[,-c(1,2)],2,max)))),
          axes = F)
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,ceiling(max(early[,c(2:ncol(early))])),100)) 
    
  boxplot(early[,-1],
          ylab="Expression UA",
          xlab="Temps autogamie",
          main=paste0("Expression early peak in ", condition),
          ylim = c(0,zoom_max),
          axes = F)
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,zoom_max,pas)) 
  
  boxplot(intermediate[,-1],
          ylab="Expression UA",
          xlab="Temps autogamie",
          main=paste0("Expression intermediate peak in ", condition),
          axes = F)
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,ceiling(max(intermediate[,c(2:ncol(intermediate))])),100)) 
  
  boxplot(intermediate[,-1],
          ylab="Expression UA",
          xlab="Temps autogamie",
          main=paste0("Expression intermediate peak in ", condition),
          ylim = c(0,zoom_max),
          axes = F)
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,zoom_max,pas)) 
  dev.off()  

}

CurveDrawer = function(tab, type, files, parameter){
    tab = na.omit(tab)
    condition = sub(paste0("_expression_table_",type,".tab"),"",files)
    temp = merge(annotation[,c(1,3)],tab,x.by ="ID",y.by = "ID")
    
    
    early = tab[temp$EXPRESSION_PROFIL=="Early peak",]
    intermediate = tab[temp$EXPRESSION_PROFIL=="Intermediate peak",]
    late = tab[temp$EXPRESSION_PROFIL=="Late peak",]
    late_i = tab[temp$EXPRESSION_PROFIL=="Late induction",]
    other = tab[which(is.element(tab$ID,setdiff(tab$ID,c(early$ID,intermediate$ID,late$ID,late_i$ID)))),]
    
    
    png(paste0("./Graph/",condition,"/Courbes_",condition,"_",type,"_",parameter,".png"),width = 800, height = 1000)
      plot(NULL, 
           xlim = c(1,(ncol(tab)-1)), ylim=c(0,ceiling(max(tab[,c(2:ncol(tab))]))),
           axes=F,
           ylab="Niveau d'expression (unite arbitraire)",
           xlab="Temps autogamie",
           main=paste0("Expression en condition ", condition))
      axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
      axis(2, at= seq(0,ceiling(max(tab[,c(2:ncol(tab))])),200))
      legend("topleft",legend=c("other","early","intermediate","late","late induct"),col=c(gray(0.5, alpha = 0.1),3,4,2,6),lwd=2)
      for (j in 1:nrow(tab)){
        lines(1:(ncol(tab)-1),late_i[j,-1],col=6,lwd=2)
        lines(1:(ncol(tab)-1),late[j,-1],col=2,lwd=2)
        lines(1:(ncol(tab)-1),intermediate[j,-1],col=4,lwd=2)
        lines(1:(ncol(tab)-1),early[j,-1],col=3,lwd=2)
        lines(1:(ncol(tab)-1),other[j,-1],col=gray(0.5, alpha = 0.1),lwd=2)
      }
    dev.off()  
    
    
    png(paste0("./Graph/",condition,"/Courbes_early-int_",condition,"_",type,"_",parameter,".png"), width = 800, height = 1000)
    plot(NULL, 
           xlim = c(1,(ncol(tab)-1)), ylim=c(0,ceiling(max(max(early[,c(2:ncol(early))],intermediate[,c(2:ncol(intermediate))])))),
           axes=F,
           ylab="Niveau d'expression (unite arbitraire)",
           xlab="Temps autogamie",
           main=paste0("Expression en condition ", condition))
      axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
      axis(2, at= seq(0,ceiling(max(max(early[,c(2:ncol(early))],intermediate[,c(2:ncol(intermediate))]))),100))
      legend("topleft",legend=c("early","intermediate"),col=c(3,4),lwd=2)
      for (j in 1:nrow(tab)){
        lines(1:(ncol(tab)-1),intermediate[j,-1],col=4,lwd=2)
        lines(1:(ncol(tab)-1),early[j,-1],col=3,lwd=2)
      }
    dev.off()  
    
  
}

CurveDEG = function(tab, type, files, parameter){
  tab = na.omit(tab)
  condition = sub(paste0("_expression_table_",type,".tab"),"",files)
  
  temp = merge(read.table("DATA/EVA_Siginificantly_DEGenes_PGM_KU80.tab", header=T, sep="\t")[,c(1,3,5)],tab,x.by ="ID",y.by = "ID")


  up = temp[temp$PGM == "UP"|temp$KU80 == "UP",-c(2,3)]
  up = up[which(is.element(up$ID,na.omit(up$ID))),]
  down = temp[temp$PGM == "DOWN"|temp$KU80 ==  "DOWN",-c(2,3)]
  down = down[which(is.element(down$ID,na.omit(down$ID))),]
  none = tab[which(is.element(tab$ID,setdiff(tab$ID,c(up$ID, down$ID)))),]

  png(paste0("./Graph/",condition,"/Courbes_deg_",condition,"_",type,"_",parameter,".png"), width = 800, height = 2000)
  par(mfrow=c(2,1))
  plot(NULL, 
       xlim = c(1,(ncol(tab)-1)), ylim=c(0,ceiling(max(tab[,c(2:ncol(tab))]))),
       axes=F,
       ylab="Niveau d'expression (unite arbitraire)",
       xlab="Temps autogamie",
       main=paste0("Expression en condition ", condition))
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,ceiling(max(tab[,c(2:ncol(tab))])),100))
  legend("topleft",legend=c("up","down","not deg"),col=c(2,4,gray(0.5, alpha = 0.1) ),lwd=2)
  for (j in 1:nrow(tab)){
    lines(1:(ncol(tab)-1),down[j,-1],col=4,lwd=2)
    lines(1:(ncol(tab)-1),up[j,-1],col=2,lwd=2)
    lines(1:(ncol(tab)-1),none[j,-1],col=gray(0.5, alpha = 0.1),lwd=2)
  }
  
  plot(NULL, 
       xlim = c(1,(ncol(tab)-1)), ylim=c(0,ceiling(max(max(up[,c(2:ncol(up))],down[,c(2:ncol(down))])))),
       axes=F,
       ylab="Niveau d'expression (unite arbitraire)",
       xlab="Temps autogamie",
       main=paste0("Expression en condition ", condition))
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,ceiling(max(max(up[,c(2:ncol(up))],down[,c(2:ncol(down))]))),100))
  legend("topleft",legend=c("up","down"),col=c(2,4),lwd=2)
  for (j in 1:nrow(tab)){
    lines(1:(ncol(tab)-1),down[j,-1],col=4,lwd=2)
    lines(1:(ncol(tab)-1),up[j,-1],col=2,lwd=2)
  }
  dev.off()  
  
  
  
}

CurveControls = function(type,  id){
  files = list.files(path = paste0("./DATA/",type))
  ctrl_c = read.table(paste0("./DATA/", type,"/",files[grep("CTIP_CTRL",files)]),h = T,sep = "\t")
  icl7 = read.table(paste0("./DATA/", type,"/",files[grep("ICL7",files)]),h = T,sep = "\t")
  nd7 = read.table(paste0("./DATA/", type,"/",files[grep("ND7",files)]),h = T,sep = "\t")
  ctrl_x = read.table(paste0("./DATA/", type,"/",files[grep("XRCC4_CTRL",files)]),h = T,sep = "\t")
  
  
  b = list(ctrl_c, icl7, nd7, ctrl_x)
  condition = c("nd7_ctip", "icl7", "nd7", "nd7_xrcc4")
  
  png(paste0("./Graph/Courbes_controles_",id,".png"),width = 1000, height = 1000 )
  par(mfrow=c(2,2))
  for (i in 1:4) {
    a = b[[i]]
    
    plot(NULL,
         xlim = c(0,(ncol(a)-1)), ylim=c(0,ceiling(max(a[a$ID==id,-1]))),
         axes=F,
         xlab = "autogamie",
         ylab = "expression",
         main=condition[i])
    axis(1,at = 0:(ncol(a)-2), labels = colnames(a)[-1])
    axis(2)
    lines(0:(ncol(a)-2),a[a$ID==id,2:ncol(a)])
  }
  dev.off()
}

Differences = function(control, rnai, lvl) {
  diff = array(data = NA, c(nrow(tab),5))
  diff = as.data.frame(diff)
  colnames(diff)=c("ID","Veg","T0","T5","T10")
  diff$ID = tab$ID
  diff$Veg=rnai$Veg-control$Veg
  diff$T0=rnai$T0-control$T0
  diff$T5=rnai$T5.5-control$T5
  diff$T10=rnai$T12.5-control$T10
  
  deg = merge(annotation, diff, by.x = "ID", by.y  = "ID")
  deg = deg[which(is.element(deg$PGM,c("UP","DOWN")) | is.element(deg$KU80,c("UP","DOWN"))),]
  
  Veg = diff$ID[abs(diff$Veg) > lvl]
  T0 = diff$ID[abs(diff$T0) > lvl]
  T5 = diff$ID[abs(diff$T5) > lvl]
  T10 = diff$ID[abs(diff$T10) > lvl]
  diff = diff[which(is.element(diff$ID,unique(c(Veg,T0,T10,T5)))),]
  
  diff = merge(annotation, diff, by.x = "ID", by.y  = "ID")
  
  a = list(deg, diff)
  names(a)=c("deg","diff")
  
  return(a)
}

CurveSimple = function(){
  png(paste0("./Graph/",condition,"/Courbes_",condition,"_",type,"_",parameter,".png"),width = 800, height = 1000)
  plot(NULL, 
       xlim = c(1,(ncol(tab)-1)), ylim=c(0,ceiling(max(tab[,c(2:ncol(tab))]))),
       axes=F,
       ylab="Niveau d'expression (unite arbitraire)",
       xlab="Temps autogamie",
       main=paste0("Expression en condition ", condition))
  axis(1,at=1:(ncol(tab)-1),labels=colnames(tab)[2:ncol(tab)],las=2)
  axis(2, at= seq(0,ceiling(max(tab[,c(2:ncol(tab))])),200))
  legend("topleft",legend=c("other","early","intermediate","late","late induct"),col=c(gray(0.5, alpha = 0.1),3,4,2,6),lwd=2)
  for (j in 1:nrow(tab)){
    lines(1:(ncol(tab)-1),late_i[j,-1],col=6,lwd=2)
    lines(1:(ncol(tab)-1),late[j,-1],col=2,lwd=2)
    lines(1:(ncol(tab)-1),intermediate[j,-1],col=4,lwd=2)
    lines(1:(ncol(tab)-1),early[j,-1],col=3,lwd=2)
    lines(1:(ncol(tab)-1),other[j,-1],col=gray(0.5, alpha = 0.1),lwd=2)
  }
  dev.off() 
}
