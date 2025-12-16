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
