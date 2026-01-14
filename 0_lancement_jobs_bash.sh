#!/bin/bash
##############################################
# Pipeline analyse méthylation / génomique
# Auteur : Besnard Florian
# Date   : 2025-16-12
##############################################

set -euo pipefail

##############################################
### PARAMÈTRES GÉNÉRAUX
##############################################

PROJECT_DIR="/espace_projets/inrae_gabi/rumigen/"
SCRIPT_DIR="${PROJECT_DIR}/BIN/BIN_FB"
DATA_DIR="${PROJECT_DIR}/DATA"
RESULTS_DIR="${PROJECT_DIR}/RES"

mkdir -p ${RESULTS_DIR}



##############################################
### 0. SCRIPT MAÎTRE
##############################################
echo "===== DÉMARRAGE PIPELINE ====="

##############################################
### 1. Traitement des données de méthylation
##############################################
echo ">>> Traitement données méthylation"
nohup Rscript350 ${SCRIPT_DIR}/01_traitement_methylation.R &

#Sortie 
#List des animaux filtrés:
# /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered
#Matrice des marques de méthylation filtrées et corrigées:
# /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor
#Matrice des marques de méthylation filtrées, corrigées et imputées par la moyenne si Na:
# /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/mat_epi_cor_imputed


##############################################
### 2. Traitement des données génomiques
##############################################

echo ">>> Traitement données génomiques"
nohup Rscript350 ${SCRIPT_DIR}/02_traitement_genomique.R /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered 66 &


##############################################
### 3. Extraction des YD
##############################################
echo ">>> Extraction des YD"
nohup Rscript350 ${SCRIPT_DIR}/03_extraction_YD.R /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered 66 &


##############################################
### 4. Création des deux échantillons apprentissage/ validation à partir d'une liste d'animaux
##############################################

echo ">>> Création des échantillons apprentissage / validation"
nohup Rscript350 ${SCRIPT_DIR}/04_creation_pop_appr_valid.R /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered 66 &


##############################################
### 4. Lancement des modèles classique et des modèles d'apprentissage / validation
##############################################

echo 

nohup Rscript350 ${SCRIPT_DIR}/05_lancement_modeles.R /espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered 66 &


##############################################
### 5. Extraction et compilation des résultats
##############################################
if [ "$RUN_RESULTS" = true ]; then
echo ">>> Compilation des résultats"
Rscript ${SCRIPT_DIR}/05_compile_results.R \
--input ${RESULTS_DIR} \
--output ${RESULTS_DIR}/summary
fi

##############################################
### FIN
##############################################
echo "===== PIPELINE TERMINÉ ====="
