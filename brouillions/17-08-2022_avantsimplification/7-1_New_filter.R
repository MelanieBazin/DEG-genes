UP_PKX = c(UP_PKX, UP_ALL = list(up_pkx))

not_UP_PKX = c(not_UP_PKX, not_up_PKX = list(not_up_pkx))

MOTIF = list(
  Motif = prom_motif$ID,
  Motif_plus = prom_motif$ID[prom_motif$STRAND == "+"],
  motif_moins = prom_motif$ID[prom_motif$STRAND == "-"])

MOTIF_uniq = list(
  Motif = unique(prom_motif$ID),
  Motif_plus = unique(prom_motif$ID[prom_motif$STRAND == "+"]),
  motif_moins = unique(prom_motif$ID[prom_motif$STRAND == "-"]))

SUPP = list(
  Inter_motif = intersect(AUTOGAMY$inter_peak,MOTIF_uniq$Motif),
  Inter_ssmotif = setdiff(AUTOGAMY$inter_peak,MOTIF_uniq$Motif),
  not_Inter_motif = setdiff(MOTIF_uniq$Motif,AUTOGAMY$inter_peak),
  Inter_UP_motif = intersect(intersect(AUTOGAMY$inter_peak,UP_PKX$UP_ALL), MOTIF_uniq$Motif),
  Inter_UP_ssmotif = setdiff(intersect(AUTOGAMY$inter_peak,UP_PKX$UP_ALL), MOTIF_uniq$Motif)
)

SUPP2 = list(
  Candidats = intersect(stdCTIP$DOWN_UP,MOTIF_uniq$Motif),
  NotCandidats = setdiff(MOTIF_uniq$Motif, stdCTIP$DOWN_UP),
  NotCandidats_inter = intersect(AUTOGAMY$inter_peak,setdiff(MOTIF_uniq$Motif, stdCTIP$DOWN_UP))
)

# Croiser les filtres
MOTIFxAUTOG = Crossinglist(MOTIF_uniq["Motif"], AUTOGAMY)

MOTIFxUP_PKX = Crossinglist(MOTIF_uniq["Motif"], UP_PKX)
MOTIFxUP_PKXxAUTOG = Crossinglist(MOTIFxUP_PKX, AUTOGAMY)

MOTIFxCTIP = Crossinglist(MOTIF_uniq["Motif"], stdCTIP)
MOTIFxCTIPxAUTOG = Crossinglist(MOTIFxCTIP, AUTOGAMY)

MOTIFxnotUP_PKX = Crossinglist(MOTIF_uniq["Motif"], not_UP_PKX)
MOTIFxnotUP_PKXxAUTOG = Crossinglist(MOTIFxnotUP_PKX, AUTOGAMY)

MOTIFxnotCTIP = Crossinglist(MOTIF_uniq["Motif"], not_stdCTIP)
MOTIFxnotCTIPxAUTOG = Crossinglist(MOTIFxnotCTIP, AUTOGAMY)