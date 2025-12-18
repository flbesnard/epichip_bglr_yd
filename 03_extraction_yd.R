#A partir d’une liste d’individu, d’un nom de projet et de la race des individus:
#/chemin/vers/Liste_individu $projet $race
#Extrait les YD de tout les individus depuis le dernier traitement de la bdir,
    =

#Sortie :
#Matrice indiv x caractères
#Liste_caractères_extrait_avec_info :
#Correspond à la liste des caractères extrait, première colonne correspondant aux headers de la matrice précédente, les autres colonnes contenant différentes informations sur le caractère.

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
})
##To be able to launch qsub:
Sys.setenv(SGE_ROOT = "/opt/sge")

# -------------------------------
# Arguments
# -------------------------------

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("Usage: Rscript 03_extraction_yd.R <file> <project> <Race>") 
}

listeANI <- args[1]
projet <- args[2]   
RACE <- args[3]
# listeANI="/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered"
# projet="test_extraction_YD"
# RACE=66

# -------------------------------
# Parameters
# -------------------------------

