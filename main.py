import time
import argparse
import numpy as np
from osgeo import gdal
import progdyn


def arg_parser():
    """ Extraction des arguments de la ligne de commande """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--opi1", required=True,
                        help="input opi 1 (tif)", type=str)
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

    verboseprint("Import des opis et du graph...")
    opi, img = progdyn.import_opi(ARGS.opi1, ARGS.opi2)
    graph = progdyn.import_graph(ARGS.graphTif)

    # récupération des pts saisis

    verboseprint("Import des points saisis...")
    emprise = progdyn.emprise_img_gdal(img)
    size = opi.shape
    list_pts = progdyn.recup_pts_saisis(ARGS.points,
                                        ARGS.graphGeojson,
                                        ARGS.ref,
                                        emprise,
                                        size)

    # snapping pts depart & arrive

    verboseprint("Snapping des points départ et arrivé sur le graph...")

    first_pt = list_pts[0]
    point_depart = progdyn.plus_proche_voisin(first_pt, 100, graph)

    last_pt = list_pts[-1]
    point_arrive = progdyn.plus_proche_voisin(last_pt, 100, graph)

    list_pts.insert(0, point_depart)
    list_pts.append(point_arrive)

    # raccort graph initial (à corriger pr data 1)

    verboseprint("Raccords début et fin du graph...")
    raccords = progdyn.raccord_graph(graph, point_depart, point_arrive)
    raccord_start, points_raccord_start = raccords[0][0], raccords[0][1]
    raccord_end, points_raccord_end = raccords[1][0], raccords[1][1]

    # calcul du meilleur cheminement tronçon par tronçon

    resol = img.GetGeoTransform()[1]
    marge = int(ARGS.marge/resol)
    lambda1 = ARGS.lambda1
    lambda2 = 1 - lambda1
    tension = ARGS.tension
    cmin = ARGS.cmin

    # # first raccord
    # print("raccord début :", points_raccord_start[0], point_depart)
    # chemin, p, masque = raccord_start, points_raccord_start, np.zeros((size))

    verboseprint("Cheminement C1 :", list_pts[0], list_pts[1])
    chemin, p, masque = progdyn.calc_cheminement(opi, list_pts[0], list_pts[1],
                                                 marge, lambda1, lambda2, tension, cmin)
    masque = np.where(masque == 1, 255., 0.)
    for k in range(1, len(list_pts)):
        p_before = p
        verboseprint(f"Cheminement C{k} :", list_pts[k-1], list_pts[k])
        c, p, m = progdyn.calc_cheminement(opi, list_pts[k-1], list_pts[k],
                                           marge, lambda1, lambda2, tension, cmin)
        chemin += c
        # nettoyage
        chemin = progdyn.nettoyage(chemin, p, p_before)
        m = np.where(m == 1, 255., 0.)
        masque += m

    # # last raccord
    # print("raccord fin :", point_arrive, points_raccord_end[-1])
    # chemin += raccord_end
    # chemin = progdyn.nettoyage(chemin, points_raccord_end, p) # nettoyage

    # graphe final (en cours de dev)

    verboseprint("Calcul du graph final...")
    graph_final = progdyn.propag_opi(chemin, graph)

    # export

    verboseprint("Export...")
    outputpath = ARGS.outputpath
    progdyn.export(chemin, outputpath+"chemin_global.tif",
                   img.GetGeoTransform(), img.GetProjection(), gdal.GDT_Byte)
    progdyn.export(masque, outputpath+"masque_global.tif",
                   img.GetGeoTransform(), img.GetProjection(), gdal.GDT_Byte)
    progdyn.export(graph_final, outputpath+"graph_final.tif",
                   img.GetGeoTransform(), img.GetProjection(), gdal.GDT_Byte)

    end = time.perf_counter()
    verboseprint(f"Finished in {end - start}s")
