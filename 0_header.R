source("../headers.R")

species="tetraurelia"
genome="ptetraurelia_mac_51"
annotation="ptetraurelia_mac_51_annotation_v2.0"
ref_pref="pt_51"
mapper="TOPHAT"

count_dir=paste0("/data/PARAMECIUM/COVERAGE/",species,"/",genome,"/",annotation,"/mRNA/htseq/")
base_mapping_dir=paste0("/data/PARAMECIUM/MAPPING/",species,"/",genome,"/RNAseq/")
#lancerà patir de P:/Equipes/ANGES/coeus
seqlength=read.table(paste0("data/PARAMECIUM/GENOMIC/",species,"/macronucleus/ANNOTATION/",genome,"/v2/",annotation,".transcript.fa.seqlength"),h=T)
rownames(seqlength)=sub("PTET.51.1.T","PTET.51.1.G",seqlength$ID)


# adapt these paramaters to your experiments
veg_color = "deepskyblue"
early_color = "darkgreen"
inter_color = "darkorange1"
late_color = "darkred"



xrcc4_time_points=c(rep("XRCC4_VEG",1),rep("XRCC4_EARLY",1),rep("XRCC4_INTER",1),rep("XRCC4_LATE",2))
xrcc4_time_point_colors=c(rep(veg_color,1),rep(early_color,1),rep(inter_color,1),rep(late_color,2))

nd7_n2_time_points=c(rep("ND7_VEG",1),rep("ND7_EARLY",1),rep("ND7_INTER",1),rep("ND7_LATE",2))
nd7_n2_time_point_colors=c(rep(veg_color,1),rep(early_color,1),rep(inter_color,1),rep(late_color,2))

nd7_n1_time_points=c(rep("ND7_VEG",1),rep("ND7_EARLY",2),rep("ND7_LATE",2))
nd7_n1_time_point_colors=c(rep(veg_color,1),rep(early_color,2),rep(late_color,2))



# config list 
xrcc4 = 
  xrcc4_labels = 
  nd7_n2 =
  
  
  config=list(
    
    XRCC4=list(DIR="XRCC4",PREFIXES=xrcc4,LABELS=xrcc4_labels,TIME_POINTs=xrcc4_time_points, TIME_COLORS=xrcc4_time_point_colors) ,
    ND7=list(DIR="ND7",PREFIXES=nd7_n2,LABELS=nd7_n2_labels,TIME_POINTs=nd7_n2_time_points, TIME_COLORS=nd7_n2_time_point_colors) 
    #,ND7_n1=list(DIR="ND7",PREFIXES=nd7_n1,LABELS=nd7_n1_labels,TIME_POINTs=nd7_n1_time_points, TIME_COLORS=nd7_n1_time_point_colors) 
    
  )


comparisons=list(
  "VEG"=c("XRCC4_VEG","ND7_VEG"),
  "EARLY"=c("XRCC4_EARLY","ND7_EARLY"),
  "INTER"=c("XRCC4_INTER","ND7_INTER"),
  "LATE"=c("XRCC4_LATE","ND7_LATE")
)
