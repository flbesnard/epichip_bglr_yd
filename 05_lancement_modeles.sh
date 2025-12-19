#!/bin/bash
set -e

# Modèles classiques :
# Modèle 1: Matrice similarité épigénomique + matrice similarité génomique 
# Modèle 2: Matrice similarité génomique seule
# Modèle 3: Matrice similarité génomique et matrice des marques de méthylation
# Modèle 1bis: Matrice similarité épigénomique + matrice similarité génomique avec apprentissage
# Modèle 2bis: Matrice similarité génomique seule avec YD de population d'apprentissage

##Script pour lancer tous les différents modèles:

trait="$1"
cd /espace_projets/inrae_gabi/rumigen/DATA/bglr/data

SCRIPT_DIR="/espace_projets/inrae_gabi/rumigen/DATA/bglr/script"
RES_DIR="/travail/fbesnard/bglr2/${trait}"

mkdir -p "$SCRIPT_DIR"
mkdir -p "$RES_DIR"


#########################################################
# 1. MODEL 1: GRM + MRM
#########################################################
cat <<EOF > ${SCRIPT_DIR}/R_script_mod1_${trait}.R
version
library(data.table)
library(BGLR)
library(dplyr)

setwd("${RES_DIR}")

ID_list=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered",header=FALSE)
mat_epi=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor_imputed")

row_ids <- mat_epi[[1]]
df_epi <- as.data.frame(mat_epi[, -1])
rownames(df_epi) <- row_ids
mat2 = as.matrix(df_epi[ID_list\$V1, ])
valid_cols <- names(which(apply(is.na(mat2),2,sum)==0))
mrm = tcrossprod(scale(mat2[,valid_cols],center=FALSE,scale=TRUE))/length(valid_cols)

grm= fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/grm")
setDF(grm)
rownames(grm)=names(grm)
grm=as.matrix(grm)

DATA_YD=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Mat_YD_r66")
DATA_YD=DATA_YD%>%filter( Id %in% ID_list\$V1)
setDF(DATA_YD)
rownames(DATA_YD)=DATA_YD\$Id

common_ids <- Reduce(intersect, list(
  ID_list\$V1,
  rownames(mrm),
  rownames(grm),
  DATA_YD\$Id
))

grm2 <- grm[common_ids, common_ids]

y_data <- DATA_YD[common_ids, "${trait}"]
names(y_data) <- common_ids

mrm2 <- mrm[common_ids, common_ids]


model1 = BGLR(
  y= y_data,
  ETA=list(
    gebv=list(K=grm2, model="RKHS"),
    met=list(K=mrm2, model="RKHS")
  ),
  nIter=100000, burnIn=25000, thin=5
)

save(model1, file="${RES_DIR}/model1.RData")
EOF

#########################################################
# 2. MODEL 2: GRM ONLY
#########################################################
cat <<EOF > ${SCRIPT_DIR}/R_script_mod2_${trait}.R
version
library(data.table)
library(BGLR)
library(dplyr)

setwd("${RES_DIR}")

ID_list=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered",header=FALSE)

grm= fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/grm")
setDF(grm)
rownames(grm)=names(grm)
grm=as.matrix(grm)

DATA_YD=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Mat_YD_r66")
DATA_YD=DATA_YD%>%filter( Id %in% ID_list\$V1)
setDF(DATA_YD)
rownames(DATA_YD)=DATA_YD\$Id


common_ids <- Reduce(intersect, list(
  ID_list\$V1,
  rownames(grm),
  DATA_YD\$Id
))

grm2 <- grm[common_ids, common_ids]

y_data <- DATA_YD[common_ids, "${trait}"]
names(y_data) <- common_ids



model2 = BGLR(
  y=y_data,
  ETA=list(gebv2=list(K=grm2, model="RKHS")),
  nIter=100000, burnIn=25000, thin=5
)

save(model2, file="${RES_DIR}/model2.RData")
EOF

#########################################################
# 3. MODEL 3: GRM + MCpG
#########################################################
cat <<EOF > ${SCRIPT_DIR}/R_script_mod3_${trait}.R
version
library(data.table)
library(BGLR)
library(dplyr)

setwd("${RES_DIR}")

ID_list=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered",header=FALSE)
mat_epi=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor_imputed")

row_ids <- mat_epi[[1]]
df_epi <- as.data.frame(mat_epi[, -1])
rownames(df_epi) <- row_ids
mat=as.matrix(df_epi)

grm= fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/grm")
setDF(grm)
rownames(grm)=names(grm)
grm=as.matrix(grm)

DATA_YD=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Mat_YD_r66")
DATA_YD=DATA_YD%>%filter( Id %in% ID_list\$V1)
setDF(DATA_YD)
rownames(DATA_YD)=DATA_YD\$Id


common_ids <- Reduce(intersect, list(
  ID_list\$V1,
  rownames(grm),
  rownames(mat),
  DATA_YD\$Id
))

grm2 <- grm[common_ids, common_ids]

y_data <- DATA_YD[common_ids, "${trait}"]
names(y_data) <- common_ids

mat2 <- mat[common_ids,]


model3 = BGLR(
  y= y_data,
  ETA=list(
    gebv3=list(K=grm2, model="RKHS"),
    metcpg=list(X=mat2, model="BayesC")
  ),
  nIter=100000, burnIn=25000, thin=5
)

save(model3, file="${RES_DIR}/model3.RData")
EOF


#########################################################
# 4. MODEL 1BIS: GRM + MRM WITH LEARNING SPLIT
#########################################################
cat <<EOF > ${SCRIPT_DIR}/R_script_mod1bis_${trait}.R
version
library(data.table)
library(BGLR)
library(dplyr)

setwd("${RES_DIR}")

ID_list=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered",header=FALSE)

training_ids=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Liste_animaux_apprentissage_bglr_rumigen",header=FALSE)
validation_ids=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Liste_animaux_validation_bglr_rumigen",header=FALSE)


mat_epi=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor_imputed")

row_ids <- mat_epi[[1]]
df_epi <- as.data.frame(mat_epi[, -1])
rownames(df_epi) <- row_ids
mat2 = as.matrix(df_epi[ID_list\$V1, ])
valid_cols <- names(which(apply(is.na(mat2),2,sum)==0))
mrm = tcrossprod(scale(mat2[,valid_cols],center=FALSE,scale=TRUE))/length(valid_cols)

grm= fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/grm")
setDF(grm)
rownames(grm)=names(grm)
grm=as.matrix(grm)

DATA_YD=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Mat_YD_r66")
DATA_YD=DATA_YD%>%filter( Id %in% ID_list\$V1)
setDF(DATA_YD)
rownames(DATA_YD)=DATA_YD\$Id

common_ids <- Reduce(intersect, list(
  ID_list\$V1,
  rownames(grm),
  rownames(mrm),
  DATA_YD\$Id
))

grm2 <- grm[common_ids, common_ids]

y_data <- DATA_YD[common_ids, "${trait}"]
names(y_data) <- common_ids

mrm2 <- mrm[common_ids, common_ids]


# keep only validation IDs that are in the model
val_ids <- intersect(validation_ids\$V1, names(y_data))

# set those phenotypes to NA
y_data[val_ids] <- NA


model1bis = BGLR(
  y= y_data,
  ETA=list(
    gebv=list(K=grm2, model="RKHS"),
    met=list(K=mrm2, model="RKHS")
  ),
  nIter=100000, burnIn=25000, thin=5
)



save(model1bis, file="${RES_DIR}/model1bis.RData")

EOF

#########################################################
# 5. MODEL 2BIS: GRM ONLY WITH LEARNING/VALIDATION SPLIT
#########################################################
cat <<EOF > ${SCRIPT_DIR}/R_script_mod2bis_${trait}.R
version
library(data.table)
library(BGLR)
library(dplyr)

setwd("${RES_DIR}")


ID_list=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered",header=FALSE)

training_ids=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Liste_animaux_apprentissage_bglr_rumigen",header=FALSE)
validation_ids=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Liste_animaux_validation_bglr_rumigen",header=FALSE)


grm= fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/grm")
setDF(grm)
rownames(grm)=names(grm)
grm=as.matrix(grm)

DATA_YD=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Mat_YD_r66")
DATA_YD=DATA_YD%>%filter( Id %in% ID_list\$V1)
setDF(DATA_YD)
rownames(DATA_YD)=DATA_YD\$Id

common_ids <- Reduce(intersect, list(
  ID_list\$V1,
  rownames(grm),
  DATA_YD\$Id
))

grm2 <- grm[common_ids, common_ids]

y_data <- DATA_YD[common_ids, "${trait}"]
names(y_data) <- common_ids

# keep only validation IDs that are in the model
val_ids <- intersect(validation_ids\$V1, names(y_data))

# set those phenotypes to NA
y_data[val_ids] <- NA


model2bis = BGLR(
  y= y_data,
  ETA=list(gebv2=list(K=grm2, model="RKHS")),
  nIter=100000, burnIn=25000, thin=5
)

save(model2bis, file="${RES_DIR}/model2bis.RData")

EOF



#########################################################
# 6. LAUNCH MODELS IN PARALLEL
#########################################################
# Create job array script
######################################
ARRAY_SCRIPT=${SCRIPT_DIR}/job_R_BGLR_array_${trait}.sh

cat <<EOF > ${ARRAY_SCRIPT}
/opt/R/4.2.1/bin/Rscript ${SCRIPT_DIR}/R_script_mod1_${trait}.R 
/opt/R/4.2.1/bin/Rscript ${SCRIPT_DIR}/R_script_mod2_${trait}.R 
/opt/R/4.2.1/bin/Rscript ${SCRIPT_DIR}/R_script_mod3_${trait}.R 
/opt/R/4.2.1/bin/Rscript ${SCRIPT_DIR}/R_script_mod1bis_${trait}.R 
/opt/R/4.2.1/bin/Rscript ${SCRIPT_DIR}/R_script_mod2bis_${trait}.R 
EOF

chmod +x ${ARRAY_SCRIPT}

#  Launch array job
######################################
echo "Submitting array job..."
JOBID=$(qarray -N BGLR_${trait} -cwd -o /travail/fbesnard/logs/ -e /travail/fbesnard/logs/ -sync n -l h_vmem=20G -q redhat8q  ${ARRAY_SCRIPT})
echo "Array job submitted: ${JOBID}"


exit










Rscript362 R_script_$trait
#Pour méthylation avec scaling par marque :
# grm=tcrossprod(scale(all_chr_typages,center = T,scale = T))/ncol(all_chr_typages)
mat_epi=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor")

# Extract row names (FR IDs)
row_ids <- mat_epi[[1]]

# Remove first column and convert to data.frame
df_epi <- as.data.frame(mat_epi[, -1])

# Set row names
rownames(df_epi) <- row_ids

mat2=as.matrix(df_epi[ID_list$V1,])


mrm=tcrossprod(scale(mat2[,names(which(apply(is.na(mat2),2,sum) == 0))],center = F,scale = T))/ncol(mat2)

mat2clean <- mat2[,names(which(apply(is.na(mat2),2,sum) == 0))]



setDF(Meta_Data_WP6_typ)

rownames(Meta_Data_WP6_typ)=Meta_Data_WP6_typ$ANIM

model1=BGLR(y=Meta_Data_WP6_typ[ID_list$V1,"YD_LAIT"],ETA = list(gebv=list(K=grm,model="RKHS"),met=list(K=mrm,model="RKHS")),nIter=1000,burnIn =100,thin = 5)#nIter=100000 et burnIn=25000
model2=BGLR(y=Meta_Data_WP6_typ[ID_list$V1,"YD_LAIT"],ETA = list(gebv=list(K=grm,model="RKHS"),met=list(X=mat2clean,model="BayesA")),nIter=1500,burnIn =500,thin = 5)#nIter=100000 et burnIn=25000


model3=BGLR(y=Meta_Data_WP6_typ[ID_list$V1,"YD_LAIT"],ETA = list(gebv=list(K=grm,model="RKHS")),nIter=1500,burnIn =500,thin = 5)#nIter=100000 et burnIn=25000


(model2$ETA$gebv$varU+sum(model2$ETA$met$varB))/(model2$varE + model2$ETA$gebv$varU + mean(model2$ETA$met$varB))


b_scale <- (model2$ETA$met$b-mean(model2$ETA$met$b))/sd(model2$ETA$met$b)
pval <- 2*pnorm(abs(b_scale),lower.tail=FALSE)
plot( -log10(pval),type="h")

str(model1)

ETA_gebv_varU.dat=fread("/travail/fbesnard/GRM/data/ETA_met_varU.dat")

ETA_gebv_varU.dat=ETA_gebv_varU.dat %>% mutate(index=row_number())

plot(ETA_gebv_varU.dat$index,ETA_gebv_varU.dat$V1)

plot(x=model1$ETA$gebv$u,y=model1$ETA$met$u)

model1$varE
model1$ETA$gebv$varU
model1$ETA$met$varU

model1$ETA$gebv$varU/(model1$varE + model1$ETA$gebv$varU + model1$ETA$met$varU)
model1$ETA$met$varU/(model1$varE + model1$ETA$gebv$varU + model1$ETA$met$varU)

(model1$ETA$gebv$varU+model1$ETA$met$varU)/(model1$varE + model1$ETA$gebv$varU + model1$ETA$met$varU)




model1$varE
model1$ETA$gebv$varU
model1$ETA$met$varU

model1$ETA$gebv$varU/(model1$varE + model1$ETA$gebv$varU + model1$ETA$met$varU)
model1$ETA$met$varU/(model1$varE + model1$ETA$gebv$varU + model1$ETA$met$varU)

(model3$ETA$gebv$varU)/(model3$varE + model3$ETA$gebv$varU )







# Effet génomique RKHS (moyennes postérieures)
g_hat <- model2$ETA$gebv$u

# Effet met BayesA : prédiction = X %*% b
m_hat <- mat2clean %*% model2$ETA$met$b

# Variances
varG <- var(as.numeric(g_hat))
varM <- var(as.numeric(m_hat))
varE <- mean(model2$varE)   # normalement un seul nombre si saveEffects=FALSE

# Héritabilité "génétique"
h2 <- (varG+varM) / (varG + varM + varE)

# Part épigénétique
prop_met <- varM / (varG + varM + varE)

# Part génétique
prop_gen <- varG / (varG + varM + varE)

cat("h² =", round(h2,4), "\n")
cat("Part met =", round(prop_met,4), "\n")
cat("Part gen =", round(prop_gen,4), "\n")
