#Permet de créer une population d’apprentissage et de validation pour une liste d’animaux donnée :
#Liste_animaux_$projet (1 colonne)
#Afin de minimiser les relations d’apparenté entre les deux échantillons, et de mimer le mieux possible une vrai évaluation génomique
#La logique est la suivante :
#-	On sépare en deux échantillons la population selon la date de naissance :
#o	80% pour l’échantillon d’apprentissage, les plus vieux
#o	20% pour l’échantillon de validation, les plus jeunes
#-	Si un père est présent dans deux échantillons alors on retire les animaux de l’échantillon ou le père est le moins représenté
#Sortie :
#Liste_animaux_apprentissage_$projet
#Liste_animaux_validation_$projet

#Liste des paramètres
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript 04_creation_pop_appr_valid.R /chemin/vers/Liste_animaux projet race")
}

chemin_liste_animaux <- args[1]
projet <- args[2]
race <- args[3]
#chemin_liste_animaux="/espace_projets/inrae_gabi/rumigen/DATA/Metadonnees_Consolidees_with_ue_2026-01-14.csv"
#projet="bglr_rumigen"
#race=66    

library(data.table)
# 1. Charger le pedigree
#Charger chemion vers liste animaux, écrire un fichier avec la première colonne et l'utiliser dans le awk pour extraction ensuite, qu'importe le header du fichier:

MORPHE_DATA=paste0("/travail/fbesnard/MORPHE/DATA/Mortalite_",race,".csv")

Pedigree <- as.data.table(
  system(
    paste0(
      "awk -F';' '",
      "NR==FNR {p[$1]=1; next} ",
      "($1 in p) {print $1 \";\" $2 \";\" $4 \";\" $10}",
      "' ",
      chemin_liste_animaux, " ",
      MORPHE_DATA
    ),
    intern = TRUE
  )
)

col_names=c( "anim","pere","mere","datenai")
setDT(Pedigree)
Pedigree[,c(col_names):= tstrsplit(V1, " ", fixed=TRUE)][,V1:=NULL]

# 2. Séparer en deux échantillons selon la date de naissance

Pedigree[,datenai:=as.Date(datenai,format="%Y-%m-%d")]


date_cutoff=quantile(Pedigree$datenai, 0.8, na.rm = TRUE, type = 1)

# 3. Création des deux échantillons
Pop_appr=Pedigree[datenai<=date_cutoff]$anim
Pop_valid=Pedigree[datenai>date_cutoff]$anim
# 4. Retirer les animaux de l’échantillon de validation dont le père est aussi dans l’échantillon d’apprentissage
Peres_appr=Pedigree[anim%in%Pop_appr]$pere
Pop_valid_final=Pedigree[(anim%in%Pop_valid)&!(pere%in%Peres_appr)]$anim
# 5. Ecriture des résultats
fwrite(data.table(Pop_appr),paste0("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Liste_animaux_apprentissage_",projet),row.names = F,col.names = F,quote=F)
fwrite(data.table(Pop_valid_final),paste0("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/Liste_animaux_validation_",projet),row.names = F,col.names = F,quote=F)
# 6. Vérification des tailles
cat("Taille échantillon apprentissage : ",length(Pop_appr),"\n")
cat("Taille échantillon validation : ",length(Pop_valid_final),"\n")
cat("Taille totale : ",length(Pop_appr)+length(Pop_valid_final),"\n")   
# 7. Vérification des pères communs
Peres_valid=Pedigree[anim%in%Pop_valid_final]$pere
Peres_communs=intersect(Peres_appr,Peres_valid)
cat("Nombre de pères communs entre les deux échantillons : ",length(Peres_communs),"\n")    
# 8. Vérification des animaux communs
Anim_communs=intersect(Pop_appr,Pop_valid_final)
cat("Nombre d'animaux communs entre les deux échantillons : ",length(Anim_communs),"\n")    
# fin
#---------------------------------------------
