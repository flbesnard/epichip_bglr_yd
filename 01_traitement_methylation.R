# Squelette de script R pour filtration et correction des donnees methylation (WP6)
# Auteur: Besnard Florian
# Date: 2025-11-17

.libPaths("/bao/lib_R3.5.0")

##########################  Configuration et librairies ##########################
# Installer si besoin: data.table, tidyverse, matrixStats, limma, broom
library(data.table)
library(matrixStats)
library(dplyr)
library(broom)
library(ggplot2)
library(lme4)

LOG <- function(msg){cat(sprintf("[%s] %s\n", Sys.time(), msg))}

DATA_DIR <- "/espace_projets/inrae_gabi/rumigen/DATA/"
TRAV_DIR <- "/espace_projets/inrae_gabi/rumigen/TRAV/"

# Chargement des fichiers initiaux
Beta_Value_WP6 <- fread(file.path(DATA_DIR, "Beta_Value_WP6.txt"))
Meta_Data_WP6 <- fread(file.path(DATA_DIR, "Metadonnees_Consolidees_20250112.csv"))
# incoherence <- fread(file.path(DATA_DIR, "incoherence_notations_tubes.csv"))


# Conventions / hypotheses:
# - Beta_Value_WP6: table CpG x individus ou individus x CpG. Adapter transposition si besoin.
# - Meta_Data_WP6: contient colonnes 'ID', 'age', 'has_genotype' (ou equivalent), 'tis_cell1','tis_cell2','tis_cell3', 'twin_flag', 'performance' (adapter aux noms reels)
# - incoherence contient colonne 'ID_problematiques'


##Nombre de jumeaux:
table(Meta_Data_WP6$twin,Meta_Data_WP6$TYP)

##A enlever les plaques à pb:
Plaque_pb=c("208455520015","208455520023","208455520014","208455520025")

# sample=c("208454820003",
# "208454820012",
# "208455520059",
# "208455520072")


# sample_pb=c("208455520043_R08C01","208455520025_R08C01")#,"208455520031","208455520043","208455520015","208455520023","208455520014","208455520025")

Meta_Data_WP6_typ=Meta_Data_WP6[TYP==1,][twin!=1,][!(substr(Sample_Name,1,12)%in% Plaque_pb)]
# Meta_Data_WP6_typ=Meta_Data_WP6_typ[substr(Sample_Name,1,12)%in% sample ]



print(paste0("On sélectionne:",
             nrow(Meta_Data_WP6_typ),
             " samples sur un total de ",
             ncol(Beta_Value_WP6)-1,
             " samples"
))



Beta_Value_WP6=Beta_Value_WP6[, c("Marque", Meta_Data_WP6_typ$Sample_Name), with = FALSE]

# 2. Enlever animaux a probleme et marques
# 2a. Marques avec trop de NA 
Beta_Value_WP6[, prop_na_mark := rowMeans(is.na(.SD)),
               .SDcols = setdiff(names(Beta_Value_WP6), "Marque")]

marks_bad_formula <- Beta_Value_WP6[prop_na_mark > 0.2, Marque]


# 2b. Individus avec trop de NA
cols_ind <- setdiff(names(Beta_Value_WP6), "Marque")
mat <- as.matrix(Beta_Value_WP6[, ..cols_ind])

prop_na_ind <- colMeans(is.na(mat))
ids_bad_formula <- names(prop_na_ind[prop_na_ind > 0.2])

# 2a. Remove BAD MARKS (rows)
Beta_Value_WP6_marks <- Beta_Value_WP6[!(Marque %in% marks_bad_formula)]

# 2b. Remove BAD INDIVIDUALS (columns)
Beta_Value_WP6_marks_ind <- Beta_Value_WP6_marks[, 
                                                 setdiff(names(Beta_Value_WP6_marks), ids_bad_formula),
                                                 with = FALSE
]

# 3. filtre pour la SD minimum:
# Extraire uniquement les valeurs numériques (en retirant la colonne Marque, et aussi prop_na_mark)
X <- subset(Beta_Value_WP6_marks_ind, select = -c(Marque))

# Calculs par marque (ligne)
mean_per_mark <- apply(X, 1, mean, na.rm = TRUE)
sd_per_mark   <- apply(X, 1, sd,   na.rm = TRUE)

df_mark <- data.frame(
  Marque = Beta_Value_WP6_marks_ind$Marque,
  Mean = mean_per_mark,
  SD = sd_per_mark
)

setDT(df_mark)
Marques_selected=df_mark[SD>0.01]$Marque

Beta_Value_WP6_marks_ind_sd=Beta_Value_WP6_marks_ind[Marque%in%Marques_selected ]


Meta_Data_WP6_filtered=Meta_Data_WP6_typ[Sample_Name %in% names(Beta_Value_WP6_marks_ind_sd)]


fwrite(Meta_Data_WP6_filtered[,.(ANIM)],"/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered",row.names = F,col.names = F,quote=F)



##Correction:
Meta_Data_WP6_filtered[,age_prel_meta:=as.numeric(as.Date(DATE_PRELEVEMENT,"%d/%m/%Y")-as.Date(date_nais,"%d/%m/%Y"))]

#Sample plate
Meta_Data_WP6_filtered[,plate_id:=as.factor(substr(Sample_Name,1,12))]
#Sample pos
Meta_Data_WP6_filtered[,plate_pos:=as.factor(substr(Sample_Name,14,19))]



META_for_correction=Meta_Data_WP6_filtered[,.(Sample_Name,ANIM,age_prel_meta,Neutro,Lympho,Mono,plate_id,plate_pos)]



# 1. Passer les beta-values au format long
DT_long = melt(
  Beta_Value_WP6_marks_ind_sd,
  id.vars = "Marque",
  variable.name = "Sample_Name",
  value.name = "Beta"
)


# 2. Merge avec les metadata
DT_long = merge(DT_long, META_for_correction,
                by = "Sample_Name",
                all.x = TRUE)


# # Sous-ensemble de DT_long avec ces marques
# DT_sample <- DT_long[Marque %in% marques_sample]
setDT(DT_long)
DT_long <- as.data.table(DT_long)

# Préparation une fois pour toutes
DT_long[, Beta_cor  := {
  m <- lmer(
    Beta ~ age_prel_meta + Neutro + Lympho + Mono +
      (1 | plate_id) + (1 | plate_pos),
    na.action = na.exclude
  )
  resid(m)
}, by = Marque]


# Transformer en matrice large
mat <- dcast(DT_long, ANIM ~ Marque, value.var = "Beta_cor", fill = NA)
write.table(DT_long,"/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/epi_cor_long",row.names=F,col.names=F,quote=F)
write.table(mat,"/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor",row.names=F,col.names=T,quote=F)


##Imputation moyenne pour NA:
mat_epi=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor")


mat_epi[1:5,1:5]

# Pour chaque colonne numérique : remplacer les NA par la moyenne de la colonne
for (j in names(mat_epi)[-1]) {   # -1 si la première colonne est un ID
  set(mat_epi,
      i = which(is.na(mat_epi[[j]])),
      j = j,
      value = mean(mat_epi[[j]], na.rm = TRUE))
}

sum(is.na(mat_epi))   # doit afficher 0
#ok

write.table(mat_epi,"/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor_imputed",row.names=F,col.names=T,quote=F)

