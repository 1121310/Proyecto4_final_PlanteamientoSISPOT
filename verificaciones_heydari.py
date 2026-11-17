# -*- coding: utf-8 -*-
"""
Verificaciones del Reto 2 ANTES de comparar con los resultados de Heydari.

A. Datos: coherencia interna de la Tabla 4 y de los datos digitalizados/transcritos.
B. Modelo: regresion contra el codigo ya validado del Reto 1 (Garces) y dualidad
   primal = dual con eolica, HVDC y racionamiento.
C. PEM: los puntos reproducen los momentos de cada bloque.
D. Monte Carlo + k-means: Beta con los momentos de cada bloque, horas, probabilidades y medias.
Uso:  python verificaciones_heydari.py
"""
import csv
import math
import os
import sys

import datos_garver as dg
import datos_heydari as dh
from tnep_heydari import (SistemaH, LineaAC, Termica, DemandaH, Caso, costo_plan, parse_plan,
                          evaluar_h, dual_lp_h, resolver_h)
from incertidumbre import segmentacion, pem_casos, pem_puntos, momentos, montecarlo_kmeans, beta_momentos

ok_total = True


def chk(nombre, cond, detalle=""):
    global ok_total
    ok_total &= bool(cond)
    print(("  OK  " if cond else " FALLA"), nombre, detalle)


# ---------------------------------------------------------------- A. datos
print("A. Datos")
chk("bloques suman 8760 h", sum(dh.N_H) == 8760)
chk("CRF(10 %, 25 a) = 0.110168", abs(dh.CRF - 0.1101681) < 1e-6, f"({dh.CRF:.6f})")
fac = 3861.60 / 20.0
chk("costos HVAC (UTP) = Tabla I de Garces / 193.08 (15 corredores)",
    all(abs(c / fac - dh.COSTO_AC[(o, r)]) < 1e-6 for (o, r, _, _, c) in dg._TABLA_I))
# Snapshot de la pagina GP-UTP (fuente_utp_garver_hvac_hvdc.csv, copiado del HTML sin retocar)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fuente_utp_garver_hvac_hvdc.csv"),
          encoding="utf-8") as _f:
    _utp = list(csv.DictReader(ln for ln in _f if not ln.startswith("#")))
_utp_ac = {(int(f["i"]), int(f["j"])): f for f in _utp if f["tipo"] == "AC"}
_utp_dc = {(int(f["i"]), int(f["j"])): f for f in _utp if f["tipo"] == "DC"}
_tab = {(o, r): (x, fl) for (o, r, x, fl, _) in dg._TABLA_I}
chk("pagina UTP: costos HVAC, reactancias, limites de flujo y circuitos existentes = codigo (15 corredores)",
    len(_utp_ac) == 15 and all(
        abs(float(_utp_ac[k]["costo_Mdolar"]) - dh.COSTO_AC[k]) < 1e-9
        and abs(float(_utp_ac[k]["x_pu"]) - _tab[k][0]) < 1e-9
        and abs(float(_utp_ac[k]["limite_flujo_MW"]) - _tab[k][1]) < 1e-9
        and (int(_utp_ac[k]["circuitos_existentes"]) == 1) == (k in dg._EXISTENTES) for k in _utp_ac))
chk("pagina UTP: costos HVDC (15, no solo el 2-6) y limites de flujo (= HVAC del corredor) = codigo",
    len(_utp_dc) == 15 and all(
        abs(float(_utp_dc[k]["costo_Mdolar"]) - dh.COSTO_DC[k]) < 1e-9
        and abs(float(_utp_dc[k]["limite_flujo_MW"]) - _tab[k][1]) < 1e-9
        and int(_utp_dc[k]["max_circuitos"]) == 2 for k in _utp_dc))
chk("precios Tabla 3 = 1.1 x precios de Garces",
    all(abs(dh.PRECIOS_GEN[b][i] - 1.1 * dg._GENS[b][1][i]) < 1e-9 for b in dh.PRECIOS_GEN for i in range(3))
    and all(abs(dh.PRECIOS_DEM[b][i] - 1.1 * dg._DEMS[b][1][i]) < 1e-9 for b in dh.PRECIOS_DEM for i in range(3)))

S2 = dh.sistema_heydari(4)
for caso, ref in dh.PAPER_TABLA_4.items():
    pa, pd = parse_plan(ref["ac"]), parse_plan(ref["dc"])
    tot, ca, cd = costo_plan(S2, pa, pd)
    chk(f"caso {caso}: costos del plan publicado = Tabla 4 ({ref['inv']} = {ref['inv_ac']} + {ref['inv_dc']})",
        abs(tot - ref["inv"]) < 1e-6 and abs(ca - ref["inv_ac"]) < 1e-6 and abs(cd - ref["inv_dc"]) < 1e-6,
        f"-> {tot:.1f} = {ca:.1f} + {cd:.1f}")
    obj = ref["esw"] - dh.CRF * ref["inv"] + (dh.BETA * ref["cvar"] if ref["cvar"] else 0.0)
    chk(f"caso {caso}: objetivo publicado = ESW - CRF*inv" + (" + 5*CVaR" if ref["cvar"] else ""),
        abs(obj - ref["obj"]) < 0.03, f"({obj:.3f} vs {ref['obj']})")

pot_eol = 500.0 + 500.0 * 1.1
e_eol = sum(nh * sum(w) / len(w) for nh, w in zip(dh.N_H, dh.VIENTO)) * pot_eol / 1e3
chk("energia eolica (datos digitalizados x 1050 MW) ~ 2213.81 GWh publicados",
    abs(e_eol / 2213.81 - 1) < 0.015, f"({e_eol:.1f} GWh, {100 * (e_eol / 2213.81 - 1):+.2f} %)")
dem_pico = 2.2 * sum(sum(dg._DEMS[b][0]) for b in dg._DEMS)
e_dem = dem_pico * sum(nh * c for nh, c in zip(dh.N_H, dh.CARGA)) / 1e3
chk("caso 1: consumo publicado entre 90 % y 100 % de la demanda de la curva de duracion",
    0.9 * e_dem <= 9627.92 <= e_dem, f"({0.9 * e_dem:.0f} <= 9627.92 <= {e_dem:.0f} GWh)")
gen_cap = sum(sum(dg._GENS[b][0]) for b in dg._GENS)
rac_est = 0.0
for t, nh in enumerate(dh.N_H):
    dmin = 0.9 * dem_pico * dh.CARGA[t]
    rac_est += nh * sum(max(0.0, dmin - gen_cap - pot_eol * w) for w in dh.VIENTO[t]) / len(dh.VIENTO[t]) / 1e3
pub = {2: 9400.21 - 7120.48 - 2213.81, 3: 9415.87 - 7138.66 - 2213.81, 4: 9396.5 - 7119.39 - 2210.24}
chk("racionamiento por suficiencia de generacion (sin red) ~ publicado en casos 2-4",
    all(abs(rac_est / v - 1) < 0.10 for v in pub.values()),
    f"(estimado {rac_est:.1f} GWh; publicado consumo-fosil-eolica = " +
    ", ".join(f"{v:.1f}" for v in pub.values()) + ")")

# ---------------------------------------------------------------- B. modelo
print("B. Modelo")
S0 = dg.sistema_garver()
SH = SistemaH("Garver-Garces", S0.barras, S0.slack,
              [LineaAC(L.o, L.r, L.x, L.fmax, L.costo / 1000.0, L.existente, L.corredor, L.copia)
               for L in S0.lineas], [], [Termica(g.bus, list(g.bloques)) for g in S0.gens], [],
              [DemandaH(d.bus, list(d.bloques), cu_factor=100.0) for d in S0.dems],
              cmax=S0.cmax, crf=S0.crf)
casos_g = [Caso("e1", 0, 8760.0, 1.0, 0.0)]
from sistema import Escenario
import tnep
esc1 = [Escenario("1", 1.0, {})]
plan_g = {(3, 5): 1, (2, 3): 1, (4, 6): 3}
a = tnep.evaluar_plan(S0, plan_g, esc1)[0]
b = evaluar_h(SH, plan_g, {}, casos_g)["esw"]
chk("regresion Reto 1: bienestar del plan de Garces con 1 escenario", abs(a - b) < 1e-9, f"({a:.6f} vs {b:.6f})")
r0 = tnep.resolver_bilevel(S0, esc1, metodo="oraculo")
r1 = resolver_h(SH, casos_g, metodo="oraculo")
chk("regresion Reto 1: optimo bilevel de 1 escenario (plan y objetivo)",
    r0.plan == r1.plan_ac and abs(r0.obj - r1.obj) < 1e-7, f"({r1.plan_txt()}, {r1.obj:.6f} vs {r0.obj:.6f})")
r2 = resolver_h(SH, casos_g, metodo="dualidad", gamma=1000.0)
chk("regresion: MILP de dualidad = oraculo (1 escenario, Garces)", abs(r2.obj - r1.obj) < 1e-7,
    f"({r2.obj:.6f})")

casos_s, grupos_s = segmentacion(dh.N_H, dh.CARGA, dh.VIENTO, True)
for caso in (2, 3, 4):
    ref = dh.PAPER_TABLA_4[caso]
    S = dh.sistema_heydari(caso)
    pa, pd = parse_plan(ref["ac"]), parse_plan(ref["dc"])
    worst, mphi = 0.0, 0.0
    for cs in [casos_s[0], casos_s[5], casos_s[12], casos_s[13], casos_s[18]]:
        res = evaluar_h(S, pa, pd, [cs])
        d, mp = dual_lp_h(S, pa, pd, cs, gamma=None)
        worst = max(worst, abs(res["sw_casos"][0] - d) / max(1.0, abs(d)))
        mphi = max(mphi, mp)
    chk(f"caso {caso}: LP primal = LP dual con el plan publicado (5 casos)", worst < 1e-9,
        f"(err rel {worst:.1e}; max|phi| = {mphi:.0f})")

# ---------------------------------------------------------------- C. PEM
print("C. PEM de Hong")
for esq in ("2m", "2m+1"):
    casos_p, info = pem_casos(dh.N_H, dh.CARGA, dh.VIENTO, esq)
    mom_ok = True
    for t, inf in enumerate(info):
        pts = inf["puntos"]
        mu, sd, l3, l4 = inf["momentos"]
        m1 = sum(p * x for x, p in pts)
        m2 = sum(p * (x - mu) ** 2 for x, p in pts)
        m3 = sum(p * (x - mu) ** 3 for x, p in pts) / sd ** 3
        mom_ok &= abs(m1 - mu) < 1e-9 and abs(m2 - sd ** 2) < 1e-9 and abs(m3 - l3) < 1e-9
        if esq == "2m+1":
            mom_ok &= abs(sum(p * (x - mu) ** 4 for x, p in pts) / sd ** 4 - l4) < 1e-9
        mom_ok &= abs(sum(p for _, p in pts) - 1) < 1e-12 and all(x >= 0 for x, _ in pts)
    chk(f"esquema {esq}: puntos reproducen media, varianza, asimetria"
        + (" y curtosis" if esq == "2m+1" else "") + "; pesos suman 1; viento >= 0", mom_ok)
    chk(f"esquema {esq}: horas ponderadas suman 8760", abs(sum(c.horas for c in casos_p) - 8760) < 1e-9)

# ---------------------------------------------------------------- D. Monte Carlo
print("D. Monte Carlo + k-means")
bm_ok = True
for v in dh.VIENTO:
    a, b = beta_momentos(v)
    mu, sd, _, _ = momentos(v)
    bm_ok &= abs(a / (a + b) - mu) < 1e-12 and abs(a * b / ((a + b) ** 2 * (a + b + 1)) - sd ** 2) < 1e-12
chk("Beta de cada bloque con la media y la varianza de los 6 niveles", bm_ok)
for K in (3, 6, 10, 20):
    for com in (True, False):
        cm, gm, inf = montecarlo_kmeans(dh.N_H, dh.CARGA, dh.VIENTO, 2000, K, 2026, com)
        err_med = max(abs(a - b) for a, b in zip(inf["media_reducida"], inf["media_niveles"]))
        chk(f"K={K} {'comonotono' if com else 'independiente'}: horas = 8760, probabilidades suman 1, "
            "media reducida = media muestral, |media - media de 6 niveles| < 0.01",
            abs(sum(c.horas for c in cm) - 8760) < 1e-6 and abs(sum(g.delta for g in gm) - 1) < 1e-12
            and max(abs(a - b) for a, b in zip(inf["media_reducida"], inf["media_muestral"])) < 1e-9
            and err_med < 0.01, f"(err. media {err_med:.4f})")

print("\nRESULTADO:", "TODO OK" if ok_total else "HAY FALLAS")
sys.exit(0 if ok_total else 1)
