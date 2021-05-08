

cluster = list(
  ICL7 = c(rep("EARLY",1),rep("INTER",2),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  ND7 = c(rep("EARLY",1),rep("INTER",2),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  PGM = c(rep("INTER",3),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  KU80c = c(rep("EARLY",1),rep("INTER",2),rep("LATE",2),rep("VERY_LATE",1),rep("VEG",1)),
  
  CTIP_CTRL = c(rep("EARLY",1),rep("INTER",3),rep("LATE",1),rep("VEG",1)),
  CTIP = c(rep("EARLY",1),rep("INTER",2),rep("LATE",1),rep("VEG",1)),
  XRCC4_CTRL = c(rep("EARLY",1),rep("INTER",2),rep("LATE",1),rep("VEG",1)),
  XRCC4 = c(rep("EARLY",1),rep("INTER",2),rep("LATE",1),rep("VEG",1))
)

veg_color = "darkorange1"
early_color = "deepskyblue"
inter_color = "chartreuse3"
late_color = "red"
very_late_color = "deeppink2"

cluster_color = list()
for(j in names(cluster)){
  vec = cluster[[j]]
  color = c()
  for (i in 1:length(vec)){
    if (vec[i] == "VEG"){
      color=c(color,veg_color)
    }else if (vec[i] == "EARLY"){
      color=c(color,early_color)
    }else if (vec[i] == "INTER"){
      color=c(color,inter_color)
    }else if (vec[i] == "LATE"){
      color=c(color,late_color)
    }else if (vec[i] == "VERY_LATE"){
      color=c(color,very_late_color)
    }
  }
  cluster_color[[j]] = color
}
