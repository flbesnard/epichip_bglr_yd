#A partir d’une liste d’individu, d’un nom de projet et de la race des individus:
#/chemin/vers/Liste_individu $projet $race
#Extrait les YD de tout les individus depuis le dernier traitement de la bdir,

#Sortie :
#Matrice indiv x caractères
#Liste_caractères_extrait_avec_info :
#Correspond à la liste des caractères extrait, première colonne correspondant aux headers de la matrice précédente, les autres colonnes contenant différentes informations sur le caractère.

################################################################################
# Script d'extraction des YD depuis la base de données
################################################################################
# Description:
#   À partir d'une liste d'individus, d'un nom de projet et de la race,
#   extrait les YD (Yield Deviations) de tous les individus depuis le dernier
#   traitement de la base de données.
#
# Usage:
#   Rscript extract_yd.R /chemin/vers/Liste_individu projet race
#
# Arguments:
#   - Liste_individu: Chemin vers le fichier contenant la liste des individus
#   - projet: Nom du projet
#   - race: Race des individus
#
# Sorties:
#   - Matrice individus x caractères Stockées dans /travail/[USER]/YD/[projet]/Mat_YD_r[race] 
#   - Liste_caractères_extrait_avec_info: métadonnées sur les caractères Stockées dans /travail/[USER]/YD/[projet]/Liste_caracteres_extrait_avec_info  
################################################################################

# ==============================================================================
# Chargement des bibliothèques
# ==============================================================================

library(data.table)

# ==============================================================================
# Récupération des arguments en ligne de commande
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

# Vérification du nombre d'arguments
if (length(args) != 3) {
  stop("Usage: Rscript extract_yd.R /chemin/vers/Liste_individu projet race")
}

chemin_liste_individu <- args[1]
projet <- args[2]
race <- args[3]
#chemin_liste_individu="/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/fr_id_filtered"
#projet="bglr_rumigen"
#race=66


# ==============================================================================
# Lecture de la liste des caractères à extraire
# ==============================================================================

Liste_car=fread("/g2b/fbesnard/GITHUB/HAPCAR/Choix_caracteres_ibl")
Liste_car=Liste_car[Race==RACE]

# ==============================================================================
# Requête: Extraction des YD depuis le dernier traitement
# ==============================================================================
# - Récupérer uniquement les données depuis le dernier traitement
date_bdir=system(" ls /bdir/bovins/resultats/laitier/fra/  | sort -rn | head -1", intern = TRUE)

donnees_yd=NULL
for (car in Liste_car$id_caractere_ctig){
  message("Extracting YD for character: ", car)
  
  fichier_yd=paste0("/bdir/bovins/resultats/laitier/fra/",date_bdir,"/",Liste_car[id_caractere_ctig==car]$Filepath,"fich_YD" )
  
  data_yd=as.data.table(system(paste0(" awk 'NR>1{if(FILENAME==ARGV[1]){p[$1]=1}; if(FILENAME==ARGV[2]){if(p[$1]==1){print $1,$6}}}' ",chemin_liste_individu," ",fichier_yd),intern = TRUE))
  col_names=c( "Id","YD")

  setDT(data_yd)
  data_yd[,c(col_names):= tstrsplit(V1, " ", fixed=TRUE)][,V1:=NULL]

#On remet dans un tableau les YD avec le nom de colonne correspondant au caractère
    setnames(data_yd,"YD",car)

    if(is.null(donnees_yd)){
      donnees_yd=data_yd
    } else {
      donnees_yd=merge(donnees_yd,data_yd,by="Id",all=TRUE)
    }
  rm( data_yd)
}


# ==============================================================================
# Ecriture des résultats
# ==============================================================================

#On écrit la matrice des YD dans travail /projet
#Création du dossier si n'existe pas
dir.create(paste0("/travail/fbesnard/YD/",projet),recursive = T,showWarnings = F)   
fwrite(donnees_yd,paste0("/travail/fbesnard/YD/",projet,"/Mat_YD_r",race),sep="\t",row.names = F,col.names = T,quote=F)  
#On écrit la liste des caractères extraits avec les infos associées
Liste_car=Liste_car[,.(id_caractere_ctig,Explication,Filepath,Race,English_name,ETG,varg,vare)]

fwrite(Liste_car,paste0("/travail/fbesnard/YD/",projet,"/Liste_caracteres_extrait_avec_info"),sep="\t",row.names = F,col.names = T,quote=F  )  


################################################################################
# Fin du script
################################################################################
