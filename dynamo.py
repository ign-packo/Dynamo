# -*- coding: utf-8 -*-
# author : IJeuffrard
# version : v.1 20/12/2023


import heapq
from osgeo import gdal
import rasterio
import numpy as np
from scipy import signal
import geopandas as gpd
# from skimage.segmentation import flood_fill
from skimage.measure import label



def calc_emprise(opi, start, end, marge):
    """ masque les opis selon une zone de recherche de taille déterminée par la marge """
    # la marge doit être strictement supérieure à 1 avec le cout de correl
    if marge <= 1:
        marge = 2
    l, c = opi.shape[1], opi.shape[2]
    i_min = min(start[0], end[0]) - marge
    i_min = max(i_min, 0)
    i_max = max(start[0], end[0]) + marge
    i_max = min(i_max, l)
    j_min = min(start[1], end[1]) - marge
    j_min = max(j_min, 0)
    j_max = max(start[1], end[1]) + marge
    j_max = min(j_max, c)
    masque = np.ones((l, c))
    masque[:] = np.nan
    masque[i_min:i_max, j_min:j_max] = 1
    opi_masque = np.multiply(opi, masque)
    masque = np.where(masque == 1, 255., 0.)
    return opi_masque, masque


def correl_5x5(A, B):
    """ calcule le cout de correlation entre une opi A et B """
    A = np.pad(A, 2, 'edge')
    B = np.pad(B, 2, 'edge')
    SA = np.pad(A, 2, 'edge')[0:-4, 2:-2] + np.pad(A, 2, 'edge')[1:-3, 2:-2] + np.pad(A, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(A, 2, 'edge')[3:-1, 2:-2] + np.pad(A, 2, 'edge')[4:, 2:-2]
    SB = np.pad(B, 2, 'edge')[0:-4, 2:-2] + np.pad(B, 2, 'edge')[1:-3, 2:-2] + np.pad(B, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(B, 2, 'edge')[3:-1, 2:-2] + np.pad(B, 2, 'edge')[4:, 2:-2]
    SAB = np.multiply(np.pad(A, 2, 'edge')[0:-4, 2:-2], np.pad(B, 2, 'edge')[0:-4, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[1:-3, 2:-2], np.pad(B, 2, 'edge')[1:-3, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[2:-2, 2:-2], np.pad(B, 2, 'edge')[2:-2, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[3:-1, 2:-2], np.pad(B, 2, 'edge')[3:-1, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[4:, 2:-2], np.pad(B, 2, 'edge')[4:, 2:-2])
    SA2 = np.multiply(np.pad(A, 2, 'edge')[0:-4, 2:-2], np.pad(A, 2, 'edge')[0:-4, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[1:-3, 2:-2], np.pad(A, 2, 'edge')[1:-3, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[2:-2, 2:-2], np.pad(A, 2, 'edge')[2:-2, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[3:-1, 2:-2], np.pad(A, 2, 'edge')[3:-1, 2:-2]) + \
        np.multiply(np.pad(A, 2, 'edge')[4:, 2:-2], np.pad(A, 2, 'edge')[4:, 2:-2])
    SB2 = np.multiply(np.pad(B, 2, 'edge')[0:-4, 2:-2], np.pad(B, 2, 'edge')[0:-4, 2:-2]) + \
        np.multiply(np.pad(B, 2, 'edge')[1:-3, 2:-2], np.pad(B, 2, 'edge')[1:-3, 2:-2]) + \
        np.multiply(np.pad(B, 2, 'edge')[2:-2, 2:-2], np.pad(B, 2, 'edge')[2:-2, 2:-2]) + \
        np.multiply(np.pad(B, 2, 'edge')[3:-1, 2:-2], np.pad(B, 2, 'edge')[3:-1, 2:-2]) + \
        np.multiply(np.pad(B, 2, 'edge')[4:, 2:-2], np.pad(B, 2, 'edge')[4:, 2:-2])
    SA = np.pad(SA, 2, 'edge')[2:-2, 0:-4] + np.pad(SA, 2, 'edge')[2:-2, 1:-3] + np.pad(SA, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(SA, 2, 'edge')[2:-2, 3:-1] + np.pad(SA, 2, 'edge')[2:-2, 4:]
    SB = np.pad(SB, 2, 'edge')[2:-2, 0:-4] + np.pad(SB, 2, 'edge')[2:-2, 1:-3] + np.pad(SB, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(SB, 2, 'edge')[2:-2, 3:-1] + np.pad(SB, 2, 'edge')[2:-2, 4:]
    SAB = np.pad(SAB, 2, 'edge')[2:-2, 0:-4] + np.pad(SAB, 2, 'edge')[2:-2, 1:-3] + np.pad(SAB, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(SAB, 2, 'edge')[2:-2, 3:-1] + np.pad(SAB, 2, 'edge')[2:-2, 4:]
    SA2 = np.pad(SA2, 2, 'edge')[2:-2, 0:-4] + np.pad(SA2, 2, 'edge')[2:-2, 1:-3] + np.pad(SA2, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(SA2, 2, 'edge')[2:-2, 3:-1] + np.pad(SA2, 2, 'edge')[2:-2, 4:]
    SB2 = np.pad(SB2, 2, 'edge')[2:-2, 0:-4] + np.pad(SB2, 2, 'edge')[2:-2, 1:-3] + np.pad(SB2, 2, 'edge')[2:-2, 2:-2] + \
        np.pad(SB2, 2, 'edge')[2:-2, 3:-1] + np.pad(SB2, 2, 'edge')[2:-2, 4:]

    SA = SA[2:-2, 2:-2]
    SB = SB[2:-2, 2:-2]
    SA2 = SA2[2:-2, 2:-2]
    SB2 = SB2[2:-2, 2:-2]
    SAB = SAB[2:-2, 2:-2]

    Cov = SAB-np.multiply(SA, SB)/25
    # print("cov : ", Cov)
    Var = np.sqrt(np.multiply(SA2-np.multiply(SA, SA)/25., SB2-np.multiply(SB, SB)/25.))
    # print("varA : ", SA2-np.multiply(SA, SA)/25.)
    # print("varB : ", SB2-np.multiply(SB, SB)/25.)
    np.seterr(divide='ignore', invalid='ignore')     # hide exception from division by Nan
    Coef = np.divide(Cov, Var)
    Coef[Var==0.]=0.
    # print("coef : ", Coef)
    return Coef


def calc_cout_init(opis, lambda1, lambda2, tension):
    """ calcule les coûts initiaux entre 2 opis par simple différence et correlation """
    cout_diff = np.absolute(opis[0].astype(np.float64)-opis[1].astype(np.float64))
    cout_correl = correl_5x5(opis[0].astype(np.float64), opis[1].astype(np.float64))
    cout_correl = np.round(100*(np.ones(cout_correl.shape) - cout_correl), 0)
    cout = lambda1*pow(cout_diff, tension) + lambda2*pow(cout_correl, tension)
    cout = cout / pow(10, tension)
    cout = np.floor(cout+0.5)
    cout[opis[0] == 0] = np.nan
    cout[opis[1] == 0] = np.nan
    return cout


def dijkstra(cout_init, start, end, cout_min=1./256.):
    """
    calcule les couts cumulés à partir d'une matrice de couts initiale,
    les couts cumulés sont ici imagés par une 'distance'
    """
    cout_init = np.pad(cout_init, 1, 'constant', constant_values=[np.nan, np.nan])
    start = (start[0]+1, start[1]+1)
    end = (end[0]+1, end[1]+1)

    n, m = cout_init.shape
    infini = float('inf')
    distance = np.ones((n, m))*[infini]     # Initialise toutes les distances à l'infini
    distance[start] = 0

    pq = [(0, start)]       # Utilise une file de priorité pour les nœuds non visités
    heapq.heapify(pq)

    voisins = [(-1, 0), (0, 1), (0, -1), (1, 0)]

    while pq:
        # Prend le nœud avec la plus petite distance dans la file de priorité
        (dist, u) = heapq.heappop(pq)

        # Si le nœud actuel est celui que l'on cherche, on peut sortir de la boucle
        if u == end:
            break

        # Parcours tous les voisins du nœud actuel
        for i in voisins:
            v = (i[0] + u[0], i[1] + u[1])
            if np.isnan(cout_init[v]):
                continue
            # Si le voisin a une distance actuelle plus grande que la distance du nœud actuel +
            # la distance entre le nœud actuel et le voisin, on met à jour la distance du voisin.
            new_distance = distance[u] + cout_init[v] + cout_min
            if distance[v] > new_distance:
                distance[v] = new_distance
                # Ajoute le voisin avec sa nouvelle distance dans la file de priorité
                heapq.heappush(pq, (new_distance, v))

    return distance[1:-1, 1:-1]


def retour(cc, debut, fin):
    """ fonction de retour pour trouver le meilleur chemin dans une matrice de couts cumulés """
    debut = [debut[0], debut[1]]
    fin = [fin[0], fin[1]]
    # attention on rajoute un bord de np.nan
    cc = np.pad(cc, pad_width=1, mode='constant', constant_values=[np.NaN, np.NaN])
    solution = np.zeros(cc.shape)
    l = fin[0] + 1      # donc on rajoute +1 en ligne
    c = fin[1] + 1      # et en colonne
    solution[l, c] = 255
    list_points_solution = []
    list_points_solution.append((l-1, c-1))
    cc[l, c] = np.nan
    while [l-1, c-1] != debut:
        c2 = c
        l2 = l
        l_voisins = [-1, 0, 0, 1] + np.array([l2])
        c_voisins = [0, 1, -1, 0] + np.array([c2])
        C_voisins = cc[l_voisins, c_voisins]
        k = np.nanargmin(C_voisins)
        l2 = l_voisins[k]
        c2 = c_voisins[k]
        if cc[l, c] < cc[l2, c2]:
            print("erreur dans la remontee")
            break
        else:
            l = l2
            c = c2
            solution[l, c] = 255
            cc[l, c] = np.nan
            list_points_solution.append((l-1, c-1))
    solution = solution[1:-1, 1:-1]
    list_points_solution.reverse()
    return solution, list_points_solution


def calc_cheminement(opi, start, end, marge, lambda1, lambda2, tension, coutmin):
    """ retourne le meilleur chemin entre les points start et end """
    opi_masque, masque = calc_emprise(opi, start, end, marge)
    couts_init = calc_cout_init(opi_masque, lambda1, lambda2, tension)
    couts_cumul = dijkstra(couts_init, start, end, coutmin)
    chemin, points_chemin = retour(np.copy(couts_cumul), start, end)
    return chemin, points_chemin, masque


def nettoyage_agregation(agregat, points_agregat, points_chemin):
    """ nettoie les boucles et aller-retour formés par l'agregation des chemins """
    # on cherche les points de divergence (soit des points bordés par 3 voisins à 255.)
    agregat = np.where(agregat >= 255., 255., 0.)
    points_chemin_array = np.array(points_chemin)
    somme_voisins = agregat[points_chemin_array[:, 0]-1, points_chemin_array[:, 1]] + \
                    agregat[points_chemin_array[:, 0], points_chemin_array[:, 1]-1] + \
                    agregat[points_chemin_array[:, 0]+1, points_chemin_array[:, 1]] + \
                    agregat[points_chemin_array[:, 0], points_chemin_array[:, 1]+1]
    index_divergence = np.argwhere(somme_voisins > 510)

    if index_divergence.shape[0] == 0:
        pass    # pas de nettoyage si pas de divergence
    else :
        # on supprime les faux points divergence (cas où le tronçon courant colle simplement le précédent)
        list_index_divergence = []
        for i in range(index_divergence.shape[0]):
            div = index_divergence[i][0]
            if tuple(points_chemin[div]) in points_agregat:
                list_index_divergence.append(div)

        if len(list_index_divergence) == 0:
            pass  # pas de nettoyage si plus de divergence
        else :
            index_last_divergence = list_index_divergence[-1]
            point_last_divergence = points_chemin[index_last_divergence]

            # on supprime les points avant la derniere divergence sur le tronçon courant
            # suffit pour éliminer les aller-retour
            points_afac = np.array(points_chemin[0:index_last_divergence])
            agregat[points_afac[:, 0], points_afac[:, 1]] = 0
            del points_chemin[0:index_last_divergence]

            # on supprime les points au delà de la divergence sur le chemin aggregé
            # nécéssaire pour éliminer les boucles
            index_last_divergence_in_points_agregat = points_agregat.index(point_last_divergence)
            points_afac = np.array(points_agregat[index_last_divergence_in_points_agregat+1:])
            agregat[points_afac[:, 0], points_afac[:, 1]] = 0
            del points_agregat[index_last_divergence_in_points_agregat:]

    return agregat, points_agregat + points_chemin


def nettoyage_intersection(chemin, liste_chemin, graph, ref):
    """
    nettoie le chemin de ses intersections avec le graph,
    en gardant la plus grande portion de chemin sans intersection de graph
    """
    # on cherche les pts qui intersectent le graph
    index_intersection = []
    for index_point in range(len(liste_chemin)):
        point = liste_chemin[index_point]
        if graph[point] == ref:
            index_intersection.append(index_point)
    if len(index_intersection) > 2:
        index_intersection = np.array(index_intersection)
        # on récupère les deux intersections voisines les plus éloignées l'une de l'autre
        # plus précisément : leur index dans la liste du chemin global
        # dist = l'écart d'index le plus grand entre deux intersections voisines
        dist = index_intersection[1:]-index_intersection[:-1]
        arg_distmax = np.argmax(dist)
        idxA = index_intersection[arg_distmax]
        idxB = index_intersection[arg_distmax + 1]
        # on supprime tous les points avant la 1ere intersection A et tout ceux apres la 2eme B
        points_avant = np.array(liste_chemin[:idxA])
        if points_avant.shape[0] > 0:
            chemin[points_avant[:, 0], points_avant[:, 1]] = 0
        points_apres = np.array(liste_chemin[idxB+1:])
        if points_apres.shape[0] > 0:
            chemin[points_apres[:, 0], points_apres[:, 1]] = 0
    return chemin, liste_chemin[idxA:idxB+1]


# def remplir_par_diffusion(chemin, graph, start, end, ref):
#     """
#     propage l'opi de référence jusqu'au nouveau chemin de mosaiquage
#     retourne le graph final
#     """
#     # calcul de la frontiere du graph initial du coté de l'autre opi
#     filtre = np.ones((3, 3))
#     res = signal.convolve2d(graph, filtre, mode='same', boundary='symm')
#     contours = np.where(res != 9*graph, graph, 0)
#     val_autre_opi = 2   # l'opi qui n'est pas la ref vaut par convention 2 mais test au cas ou
#     if ref == 2:
#         val_autre_opi = 1
#     # frontiere sert de matrice de cout initial pour remonter le graph initial avec dijkstra
#     # donc on met un cout nul sur toute la frontiere que l'on souhaite remonter
#     frontiere = np.where(contours == val_autre_opi, 0., 255.)
#     # récupération de la liste des graines avec dijkstra sur la portion de la
#     # frontiere entre les deux intersections (start et end) du graph et du chemin
#     cc_frontiere = dijkstra(frontiere, start, end, cout_min=1./256.)
#     _, liste_graine = retour(np.copy(cc_frontiere), start, end)
#     # remplissage pour les graines nécéssaires
#     graph_final = np.where(chemin >= 255., ref, graph)
#     for graine in liste_graine:
#         if graph_final[graine] != ref:
#             graph_final = flood_fill(graph_final, graine, ref, connectivity=1)
#     return graph_final

def remplir_par_diffusion(chemin, graph, verbose=False):
    """
    propage l'opi de référence jusqu'au nouveau chemin de mosaiquage
    retourne le graph final
    """
    L=label(np.maximum(graph, chemin))
    nb_labels=np.max(L)+1

    bord=np.zeros(graph.shape)
    bord[0,]=1
    bord[-1,]=1
    bord[:,0]=1
    bord[:,-1]=1
    bord[graph==0]=1

    L_gauche=np.pad(L, 1,'edge')[1:-1, 0:-2]
    L_droite=np.pad(L, 1,'edge')[1:-1, 2:  ]
    L_haut=np.pad(L, 1,'edge')[0:-2, 1:-1]
    L_bas=np.pad(L, 1,'edge')[2:  , 1:-1]

    G={}
    label_chemin=None
    for i in range(nb_labels):
        # la ou les valeurs dans le graphe initial pour cette zone (normalement une seule valeur, sauf pour la zone du chemin)
        labels_graph_init=np.unique(graph[L==i])
        if len(labels_graph_init) == 1:
            l=labels_graph_init[0]
            # on cherche les labels en contact (en 4c)
            voisins = np.unique(np.concatenate((L_gauche[L==i], L_droite[L==i], L_haut[L==i], L_bas[L==i])))    
            voisins = voisins[voisins!=i]
            G[i]={'bord': 1. in bord[L == i], 'voisins': voisins, 'label_init': l, 'label': None}
        else:
            if label_chemin is not None:
                print('Erreur, on avait deja un label pour le chemin')
                exit(1)
            label_chemin=i
    if verbose:
        print(G, label_chemin)

    # nouvelle image du graphe
    new_graph=np.zeros(graph.shape)
    # on commence par labéliser les zones limites
    # si la zone touche le bord (ou une zone hors du graphe initialement) elle n'est pas modifiée

    nb_a_labeliser=len(G)
    if verbose:
        print("Nombre de zones a labeliser : ", nb_a_labeliser)
    for i in range(nb_labels):
        if i != label_chemin and G[i]['bord']:
            new_graph[L==i]=graph[L==i]
            G[i]['label']=G[i]['label_init']
            if verbose:
                print("region "+str(i)+" label: "+str(G[i]['label_init']))
            nb_a_labeliser-=1

    # tant que certaines zones ne sont pas labeliser
    # on cherche les zones ayant un voisin labéliser 
    while nb_a_labeliser>0:
        for i in range(nb_labels):
            if i != label_chemin and G[i]['label'] is None:
                # on cherche un voisin avec un label
                for v in G[i]['voisins']:
                    if v != label_chemin:
                        l=G[v]['label']
                        if l is not None:
                            new_graph[L==i]=l
                            G[i]['label']=l
                            if verbose:
                                print("region "+str(i)+" label: "+str(l))
                            nb_a_labeliser-=1
                            break

    # il ne reste plus que les pixels du chemin
    # on peut leur attribut le label 1 ou 2 (on prend le voisin max pour éviter les zones à 0 (hors graphe))
    ng_gauche=np.pad(new_graph, 1,'edge')[1:-1, 0:-2]
    ng_droite=np.pad(new_graph, 1,'edge')[1:-1, 2:  ]
    ng_haut=np.pad(new_graph, 1,'edge')[0:-2, 1:-1]
    ng_bas=np.pad(new_graph, 1,'edge')[2:  , 1:-1]
    masque_chemin=(L==label_chemin)
    new_graph[masque_chemin]=np.maximum(np.maximum(ng_gauche,ng_droite),np.maximum(ng_haut,ng_bas))[masque_chemin]
    return new_graph

def construire_ortho(graph, opi_ref, opi2):
    """ Retourne l'ortho RVB finale construite à partir du graph final"""
    opi_ref = import_RVB(opi_ref)
    opi2 = import_RVB(opi2)
    ortho = np.zeros(opi_ref.shape)
    for chanel in range(3):
        ortho[chanel] = np.where(graph == 1, opi_ref[chanel], 0)
        ortho[chanel] = np.where(graph == 2, opi2[chanel], ortho[chanel])
    return ortho


# -----------------------
# preparation des donnees
# -----------------------


def import_graph(path):
    """ fonction d'import du graph initial (tif) """
    img = gdal.Open(path)
    graph = img.ReadAsArray()
    return graph


def import_opi(path1, path2):
    """ fonction d'import des 2 opis (tif) en 1 """
    gdal.UseExceptions()
    img1 = gdal.Open(path1)
    opi1 = img1.ReadAsArray()
    img2 = gdal.Open(path2)
    opi2 = img2.ReadAsArray()
    opi = np.zeros((2, opi1.shape[1], opi1.shape[2]))
    opi[0, :, :] = opi1[0, :, :]
    opi[1, :, :] = opi2[0, :, :]
    return opi, img1


def import_RVB(path):
    """ import image 3 canaux """
    with rasterio.open(path) as raster_multi:
        band1 = raster_multi.read(1)
        band2 = raster_multi.read(2)
        band3 = raster_multi.read(3)
    return np.array([band1, band2, band3])


def export(data, name, geoTransform, wkt, format):
    """ fonction d'export des images monocanal """
    driver = gdal.GetDriverByName("GTiff")
    outDs = driver.Create(name, data.shape[1], data.shape[0], 1, format)
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(data)
    outDs.SetGeoTransform(geoTransform)
    outDs.SetProjection(wkt)


def export_RVB(rasters, output_name, opi):
    """ import image 3 canaux """
    with rasterio.Env():
        with rasterio.open(
            output_name,
            "w",
            driver="GTiff",
            height=rasters.shape[1],
            width=rasters.shape[2],
            count=rasters.shape[0],
            dtype=rasterio.float32,
            crs=rasterio.open(opi).crs,
            transform=rasterio.open(opi).transform
        ) as out_file:
            out_file.write(rasters.astype(rasterio.float32))


def emprise_img_gdal(img):
    """ retourne l'emprise d'une image """
    a, b, c, d, e, f = img.GetGeoTransform()
    img_width, img_height = img.RasterXSize, img.RasterYSize
    emprise = np.array((a, d, a+b*img_width, d+f*img_height))
    return emprise


def conversion_coord(emprise, size, point):
    """ retourne les coordonnées image à partir des coordonnées terrain """
    sizepixX = (emprise[2]-emprise[0])/(size[2])
    sizepixY = (emprise[3]-emprise[1])/(size[1])
    i = (point[1]-emprise[1])/sizepixY
    j = (point[0]-emprise[0])/sizepixX
    return (int(i), int(j))


def recup_pts_saisis(path_saisie, path_graph, opi_ref, emprise, size):
    """ retourne la liste des points saisis pour la retouche à partir d'un fichier geojson"""
    saisie_file = open(path_saisie)
    saisie = gpd.read_file(saisie_file)
    multipts = saisie.extract_unique_points()
    pts = gpd.GeoDataFrame({'geometry': multipts}).explode(index_parts=True)
    graph_file = open(path_graph)
    graph = gpd.read_file(graph_file)
    if opi_ref not in graph["CLICHE"].values:
        raise ValueError("Erreur dans l'import des points saisis : la valeur --ref (opi de référence) " +
                         "n'existe pas pour la clé 'CLICHE' du fichier graph (geojson)")
    # on ne garde que les points qui ne sont pas sur l'opi de référence
    res_intersection = pts.overlay(graph, how='intersection')
    pts_retouche = res_intersection[res_intersection["CLICHE"] != opi_ref]
    pts_retouche = pts_retouche.get_coordinates(index_parts=True)
    list_pts = []
    for k in range(pts_retouche.shape[0]):
        x = pts_retouche.iloc[k]["x"]
        y = pts_retouche.iloc[k]["y"]
        p = conversion_coord(emprise, size, (x, y))
        list_pts.append(p)
    return list_pts


def plus_proche_voisin(pt, rayon, graph):
    """
    retourne le plus proche voisin sur le graph initial d'un pt
    (snapping du premier et dernier pts de retouche saisis)
    attention graph = tif avec 0 pour le nodata, 1 pour l'opi de ref, 2 pour l'autre opi
    """
    x = pt[0]
    y = pt[1]
    val = graph[pt[0], pt[1]]
    # calcul de la zone de recherche avec gestion des bords
    l, c = graph.shape[0], graph.shape[1]
    i_min = max(x-rayon, 0)
    i_max = min(x+rayon+1, l)
    j_min = max(y-rayon, 0)
    j_max = min(y+rayon+1, c)
    zone_recherche = graph[i_min:i_max, j_min:j_max]
    # voisins de valeur differente de celle du pt
    zone_recherche = np.where(zone_recherche == val, 0, zone_recherche)
    voisins_frontiere = np.argwhere(zone_recherche != 0)
    voisins_frontiere = voisins_frontiere + [i_min, j_min]
    # si pas de voisin sur la frontiere trouvé dans le rayon de recherche
    if len(voisins_frontiere) == 0:
        raise ValueError("Premier ou dernier point saisi trop loin de l'arc initial")
    dist = np.sqrt((voisins_frontiere[:, 0]-x)**2 + (voisins_frontiere[:, 1]-y)**2)
    ppv = voisins_frontiere[np.argmin(dist), :]
    return tuple(ppv)
