# Dynamo

Outil python de retouche de graph de mosaïquage par programmation dynamique

## Pré-requis

Environnement Python avec les libs :
from osgeo import gdal
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)
![NumPy](https://img.shields.io/badge/Numpy-777BB4?style=for-the-badge&logo=numpy&logoColor=white)
import heapq
import geopandas as gpd
![Pandas](https://img.shields.io/badge/Pandas-2C2D72?style=for-the-badge&logo=pandas&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-%230C55A5.svg?style=for-the-badge&logo=scipy&logoColor=%white)
![SciPy](https://img.shields.io/badge/SciPy-654FF0?style=for-the-badge&logo=SciPy&logoColor=white)
import argparse
import time

## Paramètres



## Exemples d'utilisation

### data_0

Lancer la commande suivante :
```
python3 main.py -i data_0/input/opi1.tif -ii data_0/input/opi2.tif -g data_0/input/graph.tif -j data_0/input/graph.geojson -p data_0/input/saisieV2.geojson -r opi2.grf -o data_0/output/
```

### data_1

Lancer la commande suivante :
```
python3 main.py -i data_1/input/opi1.tif -ii data_1/input/opi2.tif -g data_1/input/graph.tif -j data_1/input/graph.geojson -p data_1/input/saisie.geojson -r opi2.tif -o data_1/output/
```