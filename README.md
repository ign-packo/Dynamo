# Dynamo

Outil python de retouche de graph de mosaïquage par programmation dynamique

## Pré-requis

Environnement Python avec les libs suivantes :

![Gdal](https://img.shields.io/badge/GDAL-5CAE58.svg?style=for-the-badge&logo=GDAL&logoColor=white)
![Rasterio](https://img.shields.io/static/v1?label=&message=rasterio&color=blue)
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-654FF0?style=for-the-badge&logo=SciPy&logoColor=white)
![Scikit-image](https://img.shields.io/static/v1?label=&message=scikit-image&color=orange)
![GeoPandas](https://img.shields.io/static/v1?label=&message=GeoPandas&color=purple)
![Heapq](https://img.shields.io/static/v1?label=&message=heapq&color=darkblue)
![Argparse](https://img.shields.io/static/v1?label=&message=argparse&color=darkred)
![Time](https://img.shields.io/static/v1?label=&message=time&color=yellow)

## Paramètres

| Paramètre | Abréviation | Requis | Type | Valeur par défaut | Description |
| --- | --- | --- | --- | --- | --- |
| --opi1 | -i | oui | str | | chemin complet de l'opi de référence (tif) |
| --opi2 | -ii | oui | str | | chemin complet de la seconde opi (tif) |
| --graphTif | -g | oui | str | | chemin complet du graph initial (tif) |
| --graphGeojson | -j | oui | str | | chemin complet du graph initial (geojson) |
| --points | -p | oui | str | | chemin complet du fichier geojson de saisie des points |
| --ref | -r | oui | str | | opi de référence |
| --outputpath | -o | oui | str | | chemin du dossier de sortie |
| --marge | -m | non | int | 20 | marge (en mètre) determinant la zone de recherche de meilleur chemin |
| --lambda1 | -l | non | float | 0.95 | poids du coût de différence (le coût de correlation étant 1-lambda1) |
| --tension | -t | non | int | 2 | tension sur le coût initial |
| --cmin | -c | non | float | 0.0001 | coût min de passage d'un pixel au pixel voisin, donne la précision des coûts cumulés |
| --verbose | -v | non | float | False | si True, affiche les prints |

**⚠️ Conventions à respecter : ⚠️**
* Le graph initial en tif ne contient qu'une seule bande avec les valeurs 0 pour les zones de no data, **1 pour l'opi de référence**, 2 pour l'autre opi
* Le geojson détaille les points du polygone de retouche
* le paramètre --ref doit être égal à la valeur de l'opi de référence du paramètre 'CLICHE' dans le fichier geojson du graph (à discuter)

## Exemples d'utilisation

### data_0

Lancer la commande suivante :
``` bash
python3 main.py -i data_0/input/opi2.tif -ii data_0/input/opi1.tif -g data_0/input/graph.tif -j data_0/input/graph.geojson -p data_0/input/saisieV2.geojson -r opi2.grf -o data_0/output/ -v True
```

Résultat en 28.3s dont 27.1s pour le calcul des cheminements*

### data_1

Lancer la commande suivante :
``` bash
python3 main.py -i data_1/input/opi2.tif -ii data_1/input/opi1.tif -g data_1/input/graph.tif -j data_1/input/graph.geojson -p data_1/input/saisie.geojson -r opi2.tif -o data_1/output/ -v True
```

Résultat en 41.8s dont 40.3s pour le calcul des cheminements*

*Dell Inc. Precision 3561, mémoire : 32,0 Gio, processeur : 11th Gen Intel® Core™ i7-11800H @ 2.30GHz × 16


***
[![IGN](images/IGN_logo.png)](IGN_logo)
