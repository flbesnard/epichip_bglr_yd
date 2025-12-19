##Rustine pour intégrer les YD de l'INRAE


##Matrice YD
library(data.table)

# On vas chercher les YD déjà existant construit avec 03_extraction_yd.R
donnees_yd=fread(paste0("/travail/fbesnard/YD/bglr_rumigen/Mat_YD_r66"))
#On vas chercher la liste des caractères extraits
Liste_car=fread(paste0("/travail/fbesnard/YD/bglr_rumigen/Liste_caracteres_extrait_avec_info"))


#On récupère les YD extrait pas Rachel provenant de l'INRAE
data_UE_Rumigen=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/data_inrae/data_UE_Rumigen.csv")
data_ia_UE_Rumigen=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/data_inrae/ia_UE_Rumigen.csv")
data_mc_pin_Rumigen=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/data_inrae/mc_pin_Rumigen.csv")
data_nec_UE_Rumigen=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/data_inrae/nec_UE_Rumigen.csv")
data_prod_UE_Rumigen=fread("/espace_projets/inrae_gabi/rumigen/DATA/bglr/data/data_inrae/prod_UE_Rumigen.csv")  

##Liste des caractères
# Impossible pour moi de refaire des YD à partir de ces données