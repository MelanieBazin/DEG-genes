#
# Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
#
# Ce contenu est sous Licence CC BY-SA 4.0
# https://creativecommons.org/licenses/by-sa/4.0/
#

# Le script ci dessous présente les différentes étapes de la mise en 
# application du programme DESeq2 pour la recherche de gène différentiellement
# exprimés. Il reprend les explications initialement présentées par les
# auteurs de l'outil dans le document suivant :

#--------------------------------------------------
# Love, M., Huber, W., & Anders, S. (2014). 
# Beginner’s guide to using the DESeq2 package. 
# In bioRxiv. https://doi.org/10.1101/002832
#--------------------------------------------------

#--------------------------------------------------
# Installation/Chargement des librairies
#--------------------------------------------------

# La première étape consiste à charger les librairies nécessaires à la
# reproduction de ces analyses. Deux librairies sont nécessaires,
# celle qui comporte les fonctionnalités du programme DESeq2, et celle
# qui donne accès au jeu de données test.

# Pour plus d'informations concernant l'installation des librairies
# il est utile de consulter les pages suivantes :
# http://bioconductor.org/packages/release/bioc/html/DESeq2.html
# http://bioconductor.org/packages/release/data/experiment/html/parathyroidSE.html

# Chargement des librairies
library("DESeq2") 
library("parathyroidSE")

#--------------------------------------------------
# Restauration de la table de comptage
#--------------------------------------------------

# Une table de comptage est un tableau qui comporte en ligne les gènes, en
# colonne les différentes librairies (ou banques) séquencées et à 
# l'intersection ligne colonne les valeurs de comptage, i.e. nombre de "reads" 
# alignés sur un gène donné lors de l'analyse RNAseq d'une libraire donnée. 
# Ces valeurs de comptage sont des nombres entiers, supérieur ou égal à 0.

# Les commandes ci dessous permettent de recréer la table de comptage,
# à partir des informations importées via la librairie R "parathyroidSE".
# A noter que cette procédure est spécifique ici de la librairie R utilisée,
# elle n'est donc pas à prendre en compte dans une mise en application
# de DESeq2 sur un autre jeu de données.

# Importation des données
data("parathyroidGenesSE") 
# --> L'objet parathyroidGenesSE est maintenant disponible dans la session
#     R courante.

# La table de comptage est accessible en utilisant la commande suivante 
countdata <- as.matrix(assay(parathyroidGenesSE))

head(countdata)
# --> Les gènes sont en ligne et les librairies séquencées en colonne. 
#     Attention, ici les librairies ne sont pas nommées. Cela est fait en 
#     utilisant la commande suivante.
coldata <- as.data.frame(colData(parathyroidGenesSE))
colnames(countdata) <- coldata$run

head(countdata)
# --> Maintenant les colonnes sont correctement nommées.

#--------------------------------------------------
# Table de description de l'expérience (metadata)
#--------------------------------------------------

# En plus de la table de comptage, le programme DESeq2 a besoin de connaitre
# la liste des facteurs qui potentiellement influencent des données de 
# comptage, associées aux gènes. Ce sont les "metadata" de l'expérience.
# Ces informations ont déjà été extraites, elles sont disponible dans l'objet
# suivant
coldata
# --> La colonne "run" a été utilisée pour nommer les colonnes de la table
#     de comptage (un "run" correspond au séquençage d'une librairie ou banque).
# --> D'autres informations sont disponibles ("experiment", "patient",
#     "treatment", "time", "submission", "study", "sample"). Elles sont
#     utiles pour la réalisation d'une bonne estimation des paramètres de
#     la modélisation statistique de DESeq2.

# Bilan : Nous disposons des données nécessaires à l'utilisation de DESEq2 !
# Un objet de R de type "matrix" qui est la table de comptage et un objet R
# de type "data.frame" qui est l'ensemble des metadata disponibles (facteurs).

#--------------------------------------------------
# Création de l'objet de données au format DESeq2
#--------------------------------------------------

# La procédure décrite ci dessous est cette fois généralisable, quelque soit
# le jeu de données (à condition que les tables soient organisées comme
# les tables "countdata" et "coldata" présentées ci dessus).

# La commande suivante permet de créer l'objet de données DESeq2. A noter
# ici que nous allons limiter notre analyse différentielle à l'effet de 
# des deux facteurs : "patient" et "treatment".
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata, 
                                            colData   = coldata, 
                                            design    = ~ patient + treatment)

#--------------------------------------------------
# FACULTATIF : Regroupement des réplicats techniques
#--------------------------------------------------

# Dans le jeu de données utilisé ici, il existe des échantillons qui ont
# été séquencés plusieurs fois. C'est le cas de l'échantillon "SRS308873"
table(coldata[,"sample"])

# Affichage des autres facteurs pour cet échantillons
coldata[coldata$sample == "SRS308873",]
# --> Deux librairies sont associées au même échantillon, il s'agit de 
#     SRR479060 et SRR479061. 

# POINT DE VIGILANCE : Il est important de faire la distinction entre les 
# terminologies "librairies" (ou "run") et "échantillons". 
# Les échantillons correspondent au matériel biologique qui a été préparé 
# pour être étudié par séquençage RNAseq. Ainsi à partir d'un unique 
# échantillon, différentes librairies peuvent être produites et séquencées.

# Il existe une fonction pour regrouper les réplicats. Ici les librairies
# sont regroupées (run), en fonction des échantillons (sample).
ddsCollapsed <- collapseReplicates(ddsFullCountTable, 
                                   groupby = ddsFullCountTable$sample, 
                                   run     = ddsFullCountTable$run)

# Regardons le résultat...
as.data.frame(colData(ddsCollapsed))
# --> Une colonne supplémentaire a été ajoutée, elle est nommée 
#     "runsCollapses". Elle indique les librairies qui ont été regroupées,
#     pour des échantillons uniques.

# Dans la table de comptage associée, les valeurs de comptages (countdata)
# ont été sommées
head(counts(ddsFullCountTable)[,ddsFullCountTable$sample == "SRS308873"])
head(counts(ddsCollapsed)[,ddsCollapsed$sample == "SRS308873"])

# La suite du travail sera donc réalisée avec l'objet R "ddsCollapse", qui 
# est bien au format DESeq2. A noter que les données de comptages, ainsi 
# que les metadata (facteurs) peuvent être retrouvées à tout moment en
# utilisant les fonctions suivantes :
counts(ddsCollapsed)
colData(ddsCollapsed)

#--------------------------------------------------
# Réalisation de l'analyse différentielle avec DESeq2
#--------------------------------------------------

# ***** Etape 1 : Préparation des données, choix de la comparaison

# Ici, nous allons comparer les expression des gènes entre les conditions
# de traitement des cellules (colonne "treatment"), au temps 48h 
# (colonne "time").

# Un nouvel objet R est donc créé, avec seulement le temps "48h" 
dds <- ddsCollapsed[,ddsCollapsed$time == "48h"]

# Il est conseillé d'appliquer cette ligne de code, afin de restructurer
# l'organisation de l'objet DESeq2 dans le cas où tous les échantillons 
# associés à une catégorie de facteur aient été supprimés.
dds$time <- droplevels(dds$time)

# Il est également conseillé de définir la catégorie de la variable "treatment"
# qui doit être utilisée comme référence dans les calculs de log2 fold changes,
# c'est à dire placée au dénominateur des ratios.

# Ici, c'est la catégorie "Control" qui est choisie comme référence,
# pour la comparaison
dds$treatment <- relevel(dds$treatment, "Control")

# Vérifions que notre table de metadata est conforme aux choix réalisés
as.data.frame(colData(dds))
# --> Tout est correct ! 
#     La préparation des données est terminée (enfin !).

# ***** Etape 2 : Lancement de la chaîne d'analyse DESEq2

# Une seule ligne de code est nécessaire (nous utilisons ici les valeurs
# des paramètres par défaut)
dds <- DESeq(dds)
# --> Les étapes accomplies sont :
#     estimating size factors
#     estimating dispersions
#     gene-wise dispersion estimates
#     mean-dispersion relationship
#     final dispersion estimates
#     fitting model and testing

# Pour mieux comprendre la signification de ces étapes, 
# il est très utile de lire :
#--------------------------------------------------
# Love, M. I., Huber, W., & Anders, S. (2014). 
# Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. 
# Genome Biology, 15(12), 1–21. https://doi.org/10.1186/s13059-014-0550-8
#--------------------------------------------------

# POINT DE VIGILANCE : Pour cette démonstration la version de la library 
#                      DESeq2 est DESeq2_1.28.1. Ainsi, l'argument "betaPrior"
#                      est par défaut spécifié "FALSE". Les valeurs 
#                      log fold change qui seront calculées sont les valeurs
#                      MLE (et non MAP) --> Voir la publication Love et al. 
#                      Genome Biology (2014)
#                      Figure 2 pour plus de détails.

#--------------------------------------------------
# Expertise des résultats obtenus
#--------------------------------------------------

# ***** Information 1 : Paramètres pour tous les gènes issus de 
#                       l'analyse statistique DESeq2

res <- results(dds)
res
# --> POINT DE VIGILANCE : Nous avons demandé un analyse différentielle, en
#     fonction du facteur "treatment" et nous avons indiqué que la catégorie
#     "Control" de ce facteur devait être choisie comme référence. 
#     L'affichage suivant :
# "log2 fold change (MLE): treatment OHT vs Control"
#     montre que "Control" a bien été utilisée en référence, mais 
#     il est important de noter que c'est la catégorie "OHT" qui 
#     a été utilisée dans la comparaison. Attention, cela ne signifie pas que
#     la catégorie "DPN" a été ignorée, elle été prise en compte dans les 
#     estimations des paramètres des modèles de DESeq2. Les résultats sont
#     simplement présentés ici en comparant "OHT" vs "Control".

# Pour changer la présentation des résultats, il est possible d'utiliser
# la commande suivante :
res2 <- results(dds, contrast = c("treatment", "DPN", "Control"))
res2
# --> Cette fois les comparaisons présentées concernent "DPN" vs "Control",
#     ce qui a été demandé avec l'argument "contrast" de la fonction "results".
#     Noter ici que les calculs de DESeq2 n'ont pas été relancés, c'est
#     la présentation des résultats qui est ajustée ici.

# ***** Information 2 : Signification des paramètres obtenus pour chaque
#                       gene.

# Commande pour obtenir une description des paramètres compilés dans la table
# de résultats
mcols(res, use.names = TRUE)
# --> En quelques mots :
#     baseMean       : Valeur moyenne des comptages normalisés 
#                      (sur toutes les librairies)
#     log2FoldChange : Taille d'effet entre les conditions comparées
#     lfcSE          : Erreur standard associée à la taille d'effet (ci dessus)
#     stat           : Statistique du test d'hypothèse réalisé
#     pvalue         : Valeurs P calculée si H0 est vraie
#     padj           : Valeurs P corrigées, prise en compte de la problématique 
#                      des test multiples
#
# --> POINT DE VIGILANCE : Il est possible que des valeurs de pvalues soient 
#     notées NA (Not Available) dans le tableau. Cela signifie soit que toutes 
#     les valeurs de comptage (count) étaient égales à 0, ou bien qu'il y 
#     avait des valeurs extrêmes, considérées comme "outliers" vis à vis de
#     la modélisation DESeq2.

# ***** Information 3 : Graphiques "diagnostiques"

# Graphique du paramètre de dispersion
plotDispEsts(dds, ylim = c(1e-6, 1e1))
# --> Chaque gène est représenté par un point noir avec en abscisse la valeur moyenne
#     des comptages normalisés (baseMean, voir information 2) et en ordonné l'estimateur
#     de la dispersion calculé en utilisant uniquement les données des gènes individuels.
# --> La courbe rouge est calculée à partir des points noirs, elle permet de modéliser la 
#     dépendance entre les mesures de dispersion et les valeurs moyennes de comptage normalisées.
# --> Les points de couleur bleu correspondent aux gènes (les mêmes que ceux représentés en noir),
#     après la procédure de "shrinkage" de DESeq2. Cette procédure consiste à réévaluer les 
#     valeurs de dispersion des gènes, en utilisant le modèle de dépendance entre la dispersion
#     et la moyenne, estimé par la ligne rouge. Ainsi, les points de couleur bleu sont plus
#     "ressérés" autour de la ligne rouge (--> C'est la fameux "shrinkage").
# --> Les gènes qui sont entourés d'un cercle bleu, sont les gènes pour lesquels les comptages
#     ne sont pas cohérents, ils sont alors considérés comme "outliers" et exclus du "shrinkage".

# Graphique MA plot
plotMA(res, ylim = c(-1, 1))
# --> Chaque gène est représenté par un point avec en abscisse la valeur moyenne
#     des comptages normalisés (baseMean, voir information 2), et en ordonné la
#     valeur du log fold change (base 2) (log2FoldChange, voir information 2).

# Les gènes représentés par des points de couleur bleu sont des gènes pour lesquels
# la valur P ajustée (padj, voir information 3) est inférieure à 0.1. cela peut être
# changé via l'argument "alpha" de la fonction (ci dessous changé à 5%) 
plotMA(res, ylim = c(-1, 1), alpha = 0.05)
# --> Il est intéressant d'observer que les gènes avec les expressions les plus
#     élevées (baseMean grand), sont sélectionnés pour des valeurs de log fold change
#     plus faibles.

#--------------------------------------------------
# Extraire une liste de gènes candidats
#--------------------------------------------------

# La paramètres statistiques calculés avec la méthodologie DESeq2 permettent
# de sélectionner des gènes intéressants pour la suite des analyses. Ces
# gènes sont de bons candidats pour être différentiellement exprimés.

# A cette étape différents stratégies de sélection des gènes sont possibles et 
# dépendent de la sensibilité du chercheur qui réalise une analyse. Ainsi il 
# est possible :

# Filtrage 1 : sélectionner tous les gènes sont la p-value ajustée est inférieure
#              à un seuil donné (T)
T <- 0.05
DEgenes <- res[which(as.data.frame(res)$padj < T),]
dim(DEgenes)
# --> 51 gènes sont sélectionnés

# Filtrage 2 : sélectionner tous les gènes sont la p-value ajustée est inférieure
#              à un seuil donné (T1) et la valeur absolue du log2FoldChange 
#              supérieure à une valeur données (T2)
T1 <- 0.1
T2 <- 0.5
DEgenes <- res[intersect(which(as.data.frame(res)$padj < T1), 
                         which(abs(as.data.frame(res)$log2FoldChange) > T2)),]
dim(DEgenes)
# --> 11 gènes sont cette fois sélectionnés

# Et bien d'autres solutions, à vous de décider !

#--------------------------------------------------
# BONUS : Changer les paramètres de la fonction DESeq2
#--------------------------------------------------

# Utilisation des paramètres par défaut
dds <- DESeq(dds)

# La documentation présente la liste des paramètres possibles :
# DESeq(object, 
#      test = c("Wald", "LRT"), 
#      fitType = c("parametric", "local", "mean"), 
#      sfType = c("ratio", "poscounts", "iterate"),
#      betaPrior, full = design(object), reduced, quiet = FALSE,
#      minReplicatesForReplace = 7, modelMatrixType, useT = FALSE,
#      minmu = 0.5, parallel = FALSE, BPPARAM = bpparam())

# Trois paramètres sont généralement ajustés, en fonction des caractéristiques
# des jeux de données :
# fitType    --> "parametric" (défaut), "local" ou "mean"
# betaPrior  --> "FALSE" (défaut) ou "TRUE"
# TestType   --> "Wald" (défaut) ou "LRT"

# Même analyse que précédemment (valeur des paramètres par défaut)
dds <- DESeq(dds, fitType = "parametric", betaPrior = FALSE, test = "Wald")
plotDispEsts(dds, ylim = c(1e-6, 1e1))

# Changement de la procédure de correction de la dispersion
dds <- DESeq(dds, fitType = "mean", betaPrior = FALSE, test = "Wald")
plotDispEsts(dds, ylim = c(1e-6, 1e1))

# Changement de la procédure de correction des log fold change
dds <- DESeq(dds, fitType = "parametric", betaPrior = TRUE, test = "Wald")
plotDispEsts(dds, ylim = c(1e-6, 1e1))
plotMA(results(dds), ylim = c(-1, 1))

# Changement de la procédure de test
dds <- DESeq(dds, fitType = "parametric", betaPrior = FALSE, test = "LRT",
             reduced = ~ treatment)
plotDispEsts(dds, ylim = c(1e-6, 1e1))
plotMA(results(dds), ylim = c(-1, 1))

# La fin :)

# Informations de la session R (fonctionnelle au moment de la réalisation de ce tutorial)
sessionInfo()
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252
#[4] LC_NUMERIC=C                   LC_TIME=French_France.1252    

#attached base packages:
#  [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] biomaRt_2.44.4              parathyroidSE_1.26.0        DESeq2_1.28.1               SummarizedExperiment_1.18.2
#[5] DelayedArray_0.14.1         matrixStats_0.57.0          Biobase_2.48.0              GenomicRanges_1.40.0       
#[9] GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0        

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.5             locfit_1.5-9.4         lattice_0.20-41        prettyunits_1.1.1      assertthat_0.2.1      
#[6] digest_0.6.27          BiocFileCache_1.12.1   R6_2.5.0               RSQLite_2.2.1          httr_1.4.2            
#[11] ggplot2_3.3.2          pillar_1.4.6           zlibbioc_1.34.0        rlang_0.4.8            progress_1.2.2        
#[16] curl_4.3               rstudioapi_0.13        annotate_1.66.0        blob_1.2.1             Matrix_1.2-18         
#[21] splines_4.0.2          BiocParallel_1.22.0    geneplotter_1.66.0     stringr_1.4.0          RCurl_1.98-1.2        
#[26] bit_4.0.4              munsell_0.5.0          compiler_4.0.2         askpass_1.1            pkgconfig_2.0.3       
#[31] openssl_1.4.3          tidyselect_1.1.0       tibble_3.0.4           GenomeInfoDbData_1.2.3 XML_3.99-0.5          
#[36] withr_2.3.0            dbplyr_2.0.0           crayon_1.3.4           dplyr_1.0.2            rappdirs_0.3.1        
#[41] bitops_1.0-6           grid_4.0.2             xtable_1.8-4           gtable_0.3.0           lifecycle_0.2.0       
#[46] DBI_1.1.0              magrittr_1.5           scales_1.1.1           stringi_1.5.3          XVector_0.28.0        
#[51] genefilter_1.70.0      xml2_1.3.2             ellipsis_0.3.1         generics_0.1.0         vctrs_0.3.4           
#[56] RColorBrewer_1.1-2     tools_4.0.2            bit64_4.0.5            glue_1.4.2             purrr_0.3.4           
#[61] hms_0.5.3              survival_3.1-12        yaml_2.2.1             AnnotationDbi_1.50.3   colorspace_2.0-0      
#[66] memoise_1.1.0
