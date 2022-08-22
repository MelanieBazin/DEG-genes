options(stringsAsFactors = FALSE)

Function = rep("", nrow(GO_therm))
GO_therm$Function = Function

# Description
GO_therm$Function[which( GO_therm$Description=="")] = "unknown"
GO_therm$Function[which(grepl("Uncharacterized protein R106", GO_therm$Description))] = "unknown"
GO_therm$Function[which(grepl("Protein of unknown function", GO_therm$Description))] = "unknown"

GO_therm$Function[which(grepl("Chromo", GO_therm$Description))] = "DNA binding"
GO_therm$Function[which(grepl("HTH", GO_therm$Description))] = "DNA binding"
GO_therm$Function[which(grepl("C2H2", GO_therm$Description))] = "DNA binding"
GO_therm$Function[which(grepl("CRC", GO_therm$Description))] = "DNA binding"
GO_therm$Function[which(grepl("transcription repressor", GO_therm$Description))] = "DNA binding"
GO_therm$Function[which(grepl("GAGA binding", GO_therm$Description))] = "DNA binding"
GO_therm$Function[which(grepl("DNA-binding", GO_therm$Description))] = "DNA binding"

GO_therm$Function[which(grepl("Histone", GO_therm$Description))] = "Chromatin remodelling"
GO_therm$Function[which(grepl("DNA gyrase", GO_therm$Description))] = "Chromatin remodelling"
GO_therm$Function[which(grepl("DNA helicase", GO_therm$Description))] = "Chromatin remodelling"
GO_therm$Function[which(grepl("Armadillo", GO_therm$Description))] = "Chromatin remodelling"
GO_therm$Function[which(grepl("JmjC", GO_therm$Description))] = "Chromatin remodelling"

GO_therm$Function[which(grepl("DnaJ", GO_therm$Description))] = "Protein binding"
GO_therm$Function[which(grepl("WD repeat", GO_therm$Description))] = "Protein binding"
GO_therm$Function[which(grepl("WD40", GO_therm$Description))] = "Protein binding"
GO_therm$Function[which(grepl("SET", GO_therm$Description))] = "Protein binding"
GO_therm$Function[which(grepl("Leucine Rich repeats", GO_therm$Description))] = "Protein binding"

GO_therm$Function[which(grepl("PAXX", GO_therm$Description))] = "DNA repair"
GO_therm$Function[which(grepl("mismatch", GO_therm$Description))] = "DNA repair"
GO_therm$Function[which(grepl("DHH family", GO_therm$Description))] = "DNA repair"
GO_therm$Function[which(grepl("Replication factor C large subunit", GO_therm$Description))] = "DNA repair"

GO_therm$Function[which(grepl("Mago nashi protein", GO_therm$Description))] = "RNA binding"

GO_therm$Function[which(grepl("ribosomal", GO_therm$Description))] = "RNA process"
GO_therm$Function[which(grepl("CID domain", GO_therm$Description))] = "RNA process"
GO_therm$Function[which(grepl("tRNA", GO_therm$Description))] = "RNA process"
GO_therm$Function[which(grepl("Non-structural protein NS-S", GO_therm$Description))] = "RNA process"
GO_therm$Function[which(grepl("Trigger factor", GO_therm$Description))] = "RNA process"


GO_therm$Function[which(grepl("Coiled coil", GO_therm$Description))] = "Coiled coil domain"
GO_therm$Function[which(grepl("PiggyBac", GO_therm$Description))] = "PiggyBac element"
GO_therm$Function[which(grepl("protease", GO_therm$Description))] = "Protease"
GO_therm$Function[which(grepl("Peptidase", GO_therm$Description))] = "Protease"
GO_therm$Function[which(grepl("metalloprotease", GO_therm$Description))] = "Protease"
GO_therm$Function[which(grepl("SUMO", GO_therm$Description))] = "Ubiquitinylation"

GO_therm$Function[which(grepl("emp24/gp25L/p24 family/GOLD", GO_therm$Description))] = "Membrane protein"
GO_therm$Function[which(grepl("receptor", GO_therm$Description))] = "Membrane protein"
GO_therm$Function[which(grepl("D-arabinono-1,4-lactone oxidase", GO_therm$Description))] = "Membrane protein"
GO_therm$Function[which(grepl("Sodium channel", GO_therm$Description))] = "Membrane protein"
GO_therm$Function[which(grepl("Glycerol-3-phosphate acyltransferase", GO_therm$Description))] = "Membrane protein"

GO_therm$Function[which(grepl("P200", GO_therm$Description))] = "Transport"
GO_therm$Function[which(grepl("Exocyst complex component 6", GO_therm$Description))] = "Transport"
GO_therm$Function[which(grepl("Dynein", GO_therm$Description))] = "Transport"

# GO_Molecular.Function
GO_therm$Function[which(GO_therm$GO.Term.descriptions..Molecular.Function.=="DNA binding")] = "DNA binding"
GO_therm$Function[which(grepl("zinc ion binding", GO_therm$GO.Term.descriptions..Molecular.Function.))] = "DNA binding"
GO_therm$Function[which(grepl("nucleic acid binding", GO_therm$GO.Term.descriptions..Molecular.Function.))] = "DNA binding"

GO_therm$Function[which(grepl("chromatin binding", GO_therm$GO.Term.descriptions..Molecular.Function.))] = "Chromatin binding"

GO_therm$Function[which(grepl("RNA binding", GO_therm$GO.Term.descriptions..Molecular.Function.))] = "RNA binding"

# GO_Biological.Process
GO_therm$Function[which(grepl("repair", GO_therm$GO.Term.descriptions..Biological.Process.))] = "DNA repair"
GO_therm$Function[which(grepl("chromosome organization", GO_therm$GO.Term.descriptions..Biological.Process.))] = "Chromatin remodelling"
GO_therm$Function[which(grepl("cell cycle", GO_therm$GO.Term.descriptions..Biological.Process.))] = "Cell cycle"

# Final
GO_therm$Function[which(GO_therm$Function=="")] = "Other process"
write.table(GO_therm, paste0(path,"GO_therm/resultsGO_BioMart2_treated.txt"), sep = "\t")

write.table(table(GO_therm$Function), paste0(path,"GO_therm/resultsGO_BioMart2_table.txt"), sep = "\t")
