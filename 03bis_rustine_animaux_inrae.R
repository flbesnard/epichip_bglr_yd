##Rustine pour intégrer les YD de l'INRAE


##Matrice YD
library(data.table)
library(ggplot2)

#On vas chercher la liste des caractères extraits
Liste_car=fread(paste0("/travail/fbesnard/YD/bglr_rumigen/Liste_caracteres_extrait_avec_info"))


# Les données brutes et les programmes (sas) pour produire les fichiers indiqués sont sous /home/dboichard/perf_ue_inrae.
# Les tables extraites par Rachel de Margau sont les fichier csv (animal, production, ia). Je n’ai pas utilisé les autres pour l’instant (nb : voulez vous l’état corporel et les mammites cliniques ?)

# Le fichier production a été créé par le pgm prep_donnees_ue.sas
# Le fichier ferti a été créé par le pgm prep_fert_ue.sas

# Vous trouverez les données des unités expé INRAE dans l’espace /espace_projet/inrae_gabi/rumigen/DATA
# Le fichier res_ue_prod contient les données pour le lait, la matière grasse, la matière protéique, le tb, le tp et les cellules.
# Les unités sont des kg de lait, mg, mp, des g/kg pour le tb et le tp, et des points de scores cellulaires (sur échelle log2)

#Lecture des fichiers de données brutes
res_ue_prod=fread("/espace_projets/inrae_gabi/rumigen/DATA/res_ue_prod")
res_ue_fert=fread("/espace_projets/inrae_gabi/rumigen/DATA/res_ue_fert")    

#Replace "." par NA dans les deux datasets
res_ue_prod[res_ue_prod == "."] <- NA
res_ue_fert[res_ue_fert == "."] <- NA


# On vas chercher les YD déjà existant construit avec 03_extraction_yd.R
donnees_yd=fread(paste0("/travail/fbesnard/YD/bglr_rumigen/Mat_YD_r66"))

#Table de correspondance entre caractère UE et caractère YD:
corresp <- data.table(
  source = c(
    "res_ue_prod", "res_ue_prod", "res_ue_prod", "res_ue_prod", "res_ue_prod", "res_ue_prod",
    "res_ue_fert", "res_ue_fert"
  ),
  nom_ue = c(
    "lait", "qmg", "qmp", "tb", "tp", "cell",
    "fertv", "fertg"
  ),
  nom_yd = c(
    "LAIT", "QMG", "QMP", "TB", "TP", "CEL",
    "FERV", "FERG"
  )
)

#On vas comparer les distributions des YD existants et des nouveaux, pour chaque caractère et on regarde les moyennes et écarts types
for (i in 1:nrow(corresp)){
  nom_yd_car=corresp$nom_yd[i]
  nom_yd_ue=corresp$nom_ue[i]

  
#On vas chercher les données ue correspondante
    if (corresp$source[i]=="res_ue_prod"){
        donnees_ue=res_ue_prod[, .(anim, get(nom_yd_ue),ue,annee) ]
    } else if (corresp$source[i]=="res_ue_fert"){
        donnees_ue=res_ue_fert[, .(anim, get(nom_yd_ue),ue,annee) ]
    }
    # remove NA
    donnees_ue=donnees_ue[!is.na(V2), .(Id=anim, V2=as.numeric(V2)) ]

  #On vas chercher les données yd correspondantes
    donnees_yd_car=donnees_yd[, .(Id, get(nom_yd_car)) ]
    

    #On fait un histogramme qui montre la distribution, la moyenne et l'écart type pour les deux jeu de données
    # calculs préalables (plus propre que dans aes)
    m_ue  <- mean(donnees_ue[, V2], na.rm = TRUE)
    sd_ue <- sd(donnees_ue[, V2], na.rm = TRUE)
    
    m_yd  <- mean(donnees_yd_car[, V2], na.rm = TRUE)
    sd_yd <- sd(donnees_yd_car[, V2], na.rm = TRUE)
    print(
    ggplot() +
      geom_histogram(data = donnees_ue, aes(x = V2),
                     fill = "blue", alpha = 0.5, bins = 50) +
      geom_vline(xintercept = m_ue,
                 color = "blue", linetype = "dashed", size = 1) +
      geom_vline(xintercept = c(m_ue - sd_ue, m_ue + sd_ue),
                 color = "blue", linetype = "dotted", size = 0.8) +
      
      geom_histogram(data = donnees_yd_car, aes(x = V2),
                     fill = "red", alpha = 0.5, bins = 50) +
      geom_vline(xintercept = m_yd,
                 color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = c(m_yd - sd_yd, m_yd + sd_yd),
                 color = "red", linetype = "dotted", size = 0.8) +
      
      labs(
        title = paste0("Comparaison des distributions pour le caractère ", nom_yd_car),
        x = "Valeur",
        y = "Fréquence"
      ) +
      theme_minimal()
)
    
   ##On enlève les outliers dans la partie prod:
    
    if (corresp$source[i]=="res_ue_prod"){
     donnees_ue=donnees_ue[V2 >= (m_ue - 1.96 * sd_ue) & V2 <= (m_ue + 1.96 * sd_ue)]

      
      #On fait une transformation pour mimer la loi normale des YD
    } else if (corresp$source[i]=="res_ue_fert"){
        #Transformation logit pour la fertilité
      donnees_ue[, V2 := logit( V2 / 100 ) ]

      
    }
    
}

