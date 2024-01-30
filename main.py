# -*- coding: utf-8 -*-
# author : IJeuffrard
# version : v.1 20/12/2023


import time
import argparse
from osgeo import gdal
import dynamo as dm


def arg_parser():
    """ Extraction des arguments de la ligne de commande """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--opi1", required=True,
                        help="input opi de reference (tif)", type=str)
    parser.add_argument("-ii", "--opi2", required=True,
                        help="input opi 2 (tif)", type=str)
    parser.add_argument("-g", "--graphTif", required=True,
                        help="input graph (tif)", type=str)
    parser.add_argument("-j", "--graphGeojson", required=True,
                        help="input graph (geojson)", type=str)
    parser.add_argument("-p", "--points", required=True,
                        help="points saisis (geojson)", type=str)
    parser.add_argument("-r", "--ref", required=True,
                        help="nom de l'opi de reference dans le geojson",
                        type=str)
    parser.add_argument("-o", "--outputpath", required=True,
                        help="outputpath", type=str)

    parser.add_argument("-m", "--marge", required=False,
                        help="marge (en mètre)",
                        type=int, default=20)
    parser.add_argument("-l", "--lambda1", required=False,
                        help="poids du cout de difference",
                        type=float, default=0.95)
    parser.add_argument("-t", "--tension", required=False,
                        help="tension",
                        type=int, default=2)
    parser.add_argument("-c", "--cmin", required=False,
                        help="cout min de passage d'un pixel au pixel voisin",
                        type=float, default=0.0001)

    parser.add_argument("-v", "--verbose",
                        help="verbose", type=bool, default=False)

    return parser.parse_args()


ARGS = arg_parser()


if __name__ == "__main__":

    verboseprint = print if ARGS.verbose else lambda *a, **k: None

    verboseprint("Arguments: ", ARGS)

    start = time.perf_counter()

    # import opis

    verboseprint("* Import des opis et du graph...")
    opi, img = dm.import_opi(ARGS.opi1, ARGS.opi2)
    graph = dm.import_graph(ARGS.graphTif)

    # récupération des pts saisis

    verboseprint("* Import des points saisis...")
    emprise = dm.emprise_img_gdal(img)
    size = opi.shape
    list_pts = dm.recup_pts_saisis(ARGS.points,
                                   ARGS.graphGeojson,
                                   ARGS.ref,
                                   emprise,
                                   size)

    # snapping pts depart & arrive

    verboseprint("* Snapping des points départ et arrivé sur le graph...")

    first_pt = list_pts[0]
    point_depart = dm.plus_proche_voisin(first_pt, 100, graph)

    last_pt = list_pts[-1]
    point_arrive = dm.plus_proche_voisin(last_pt, 100, graph)

    list_pts.insert(0, point_depart)
    list_pts.append(point_arrive)

    # calcul du meilleur cheminement tronçon par tronçon

    tic = time.perf_counter()
    verboseprint("* Calcul nouveau chemin de mosaïquage...")

    resol = img.GetGeoTransform()[1]
    marge = int(ARGS.marge/resol)
    lambda1 = ARGS.lambda1
    lambda2 = 1 - lambda1
    tension = ARGS.tension
    cmin = ARGS.cmin

    verboseprint("    Cheminement C1 :", list_pts[0], list_pts[1])
    chemin, points_chemin, masque = dm.calc_cheminement(opi, list_pts[0],
                                                        list_pts[1], marge,
                                                        lambda1, lambda2,
                                                        tension, cmin)
    chemin_global = chemin
    masque_global = masque
    points_chemin_global = points_chemin
    for k in range(2, len(list_pts)):
        verboseprint(f"    Cheminement C{k} :", list_pts[k-1], list_pts[k])
        chemin, points_chemin, masque = dm.calc_cheminement(opi, list_pts[k-1],
                                                            list_pts[k], marge,
                                                            lambda1, lambda2,
                                                            tension, cmin)
        chemin_global += chemin
        masque_global += masque
        # nettoyage
        chemin_global, points_chemin_global = dm.nettoyage_agregation(chemin_global,
                                                                      points_chemin_global,
                                                                      points_chemin)

    toc = time.perf_counter()
    verboseprint(f"{toc - tic}s")

    # nettoyage avant/apres intersections graph

    tic = time.perf_counter()
    verboseprint("* Nettoyage des intersections avec le graph...")
    REF = 1     # attention convention à respecter
    chemin_global, points_chemin_global = dm.nettoyage_intersection(chemin_global,
                                                                    points_chemin_global,
                                                                    graph, REF)
    toc = time.perf_counter()
    verboseprint(f"{toc - tic}s")

    # graph final et ortho

    tic = time.perf_counter()
    verboseprint("* Calcul du graph/ortho final...")
    graph_final = dm.remplir_par_diffusion(chemin_global,
                                           graph)
    ortho = dm.construire_ortho(graph_final, ARGS.opi1, ARGS.opi2)
    toc = time.perf_counter()
    verboseprint(f"{toc - tic}s")

    # export

    verboseprint("Export...")
    outputpath = ARGS.outputpath
    dm.export(chemin_global, outputpath+"chemin_global.tif",
              img.GetGeoTransform(), img.GetProjection(), gdal.GDT_Byte)
    dm.export(masque_global, outputpath+"masque_global.tif",
              img.GetGeoTransform(), img.GetProjection(), gdal.GDT_Byte)
    dm.export(graph_final, outputpath+"graph_final.tif",
              img.GetGeoTransform(), img.GetProjection(), gdal.GDT_Byte)
    dm.export_RVB(ortho, outputpath+"ortho.tif", ARGS.opi1)

    end = time.perf_counter()
    verboseprint(f"Finished in {end - start}s")
