"""Reto 2, metodologia 3: robustez frente a la semilla del Monte Carlo.

Repite la reduccion comonotona con K = 20 para cinco semillas y, en cada una, resuelve el
plan optimo de los casos 2-4 y lo evalua con los 30 mercados y los 6 escenarios de riesgo
de referencia. Uso:  python reto2_semillas.py [caso]      (sin argumento: casos 2, 3 y 4)
Con HiGHS:  TNEP_SOLVER=highs python reto2_semillas.py
"""
import json
import os
import sys

import datos_heydari as dh
from incertidumbre import montecarlo_kmeans
from reto2_garver import preparar, MC_N
from tnep_heydari import resolver_h, evaluar_h

SEMILLAS = (2026, 1, 2, 3, 4)
K = 20
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resultados")


def correr(caso):
    S, seg, grupos, beta = preparar(caso)
    base = resolver_h(S, seg, grupos, dh.ALFA, beta, metodo="oraculo")
    base_obj = evaluar_h(S, base.plan_ac, base.plan_dc, seg, grupos, dh.ALFA, beta)["obj"]
    filas = []
    for sem in SEMILLAS:
        cm, gm, _ = montecarlo_kmeans(dh.N_H, dh.CARGA, dh.VIENTO, MC_N, K, sem, True)
        r = resolver_h(S, cm, gm, dh.ALFA, beta, metodo="oraculo")
        ref = evaluar_h(S, r.plan_ac, r.plan_dc, seg, grupos, dh.ALFA, beta)["obj"]
        filas.append(dict(semilla=sem, plan=r.plan_txt(), inv=r.inv_total, regret=base_obj - ref,
                          igual=(r.plan_ac == base.plan_ac and r.plan_dc == base.plan_dc), t=r.t))
        print(f"caso {caso} semilla {sem}: {r.plan_txt()} regret {base_obj - ref:.3f}", flush=True)
    return dict(plan_ref=base.plan_txt(), inv_ref=base.inv_total, filas=filas)


if __name__ == "__main__":
    casos = [int(sys.argv[1])] if len(sys.argv) > 1 else [2, 3, 4]
    for c in casos:
        out = correr(c)
        with open(os.path.join(RES, f"reto2_semillas_caso{c}.json"), "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1, ensure_ascii=False)
