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
### FLAGS D'EXÉCUTION
##############################################
RUN_METHYLATION=true
RUN_GENOMIQUE=true
RUN_YD=true
RUN_MODELS=true
RUN_VALIDATION=true
RUN_RESULTS=true


##############################################
### 0. SCRIPT MAÎTRE
##############################################
echo "===== DÉMARRAGE PIPELINE ====="

##############################################
### 1. Traitement des données de méthylation
##############################################
if [ "$RUN_METHYLATION" = true ]; then
echo ">>> Traitement données méthylation"
Rscript ${SCRIPT_DIR}/01_traitement_methylation.R \
--input ${DATA_DIR}/methylation \
--output ${DATA_DIR}/methylation_processed
fi

##############################################
### 2. Traitement des données génomiques
##############################################
if [ "$RUN_GENOMIQUE" = true ]; then
echo ">>> Traitement données génomiques"
Rscript ${SCRIPT_DIR}/02_traitement_genomique.R \
--input ${DATA_DIR}/genomique \
--output ${DATA_DIR}/genomique_processed
fi

##############################################
### 3. Extraction des YD
##############################################
if [ "$RUN_YD" = true ]; then
echo ">>> Extraction des YD"
Rscript ${SCRIPT_DIR}/03_extraction_YD.R \
--phenotypes ${DATA_DIR}/phenotypes \
--output ${DATA_DIR}/YD
fi

##############################################
### 4. Lancement des modèles
##############################################
if [ "$RUN_MODELS" = true ]; then
echo ">>> Lancement modèles classiques"

# Modèle 1 : ERM + GRM
bash ${SCRIPT_DIR}/models/run_model1_ERM_GRM.sh

# Modèle 2 : GRM seule
bash ${SCRIPT_DIR}/models/run_model2_GRM.sh

# Modèle 3 : GRM + marques de méthylation
bash ${SCRIPT_DIR}/models/run_model3_GRM_METH.sh
fi

##############################################
### 4 bis. Modèles apprentissage / validation
##############################################
if [ "$RUN_VALIDATION" = true ]; then
echo ">>> Lancement modèles apprentissage / validation"

# Modèle 1 bis : ERM + GRM avec apprentissage
bash ${SCRIPT_DIR}/models/run_model1bis_learning.sh

# Modèle 2 bis : GRM seule avec YD apprentissage
bash ${SCRIPT_DIR}/models/run_model2bis_learning.sh
fi

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
