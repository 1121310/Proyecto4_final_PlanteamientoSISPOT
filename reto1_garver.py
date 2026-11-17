# -*- coding: utf-8 -*-
"""
RETO 1 - Validacion de los resultados de Garces et al. (2009) en el sistema de
Garver de 6 barras: Tabla IV (bilevel vs minimizacion de costos), fila
'basica y centralizada' de la tesis, y Tabla V (contingencias simples).

Uso:  python reto1_garver.py
Salida: resultados/garver_resultados.md  y  resultados/garver_resultados.json
"""
import json
import os
import time

import datos_garver as G
import tnep
from sistema import parse_plan, plan_a_texto

DIR_RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resultados")
os.makedirs(DIR_RES, exist_ok=True)
MD = []


def p(texto=""):
    print(texto)
    MD.append(texto)


def mismo_plan(a, b):
    return sorted(a.items()) == sorted(b.items())


def main():
    S = G.sistema_garver()
    esc = G.escenarios_base()
    out = {}
    p("# Reto 1 - Sistema de Garver: reproduccion de Garces et al. (2009)\n")
    p(f"Candidatos binarios: {len(S.candidatas)}; escenarios: {len(esc)}; sigma = {S.sigma:.0f} h; "
      f"crf (objetivo) = {S.crf}; presupuesto = {S.cmax} M EUR; M = Gamma = {S.M:.0f}; Sbase = {S.Sbase:.0f} MVA.\n")

    # ------------------------------------------------------------------
    p("## 1. Modelo bilevel (Tabla IV, primera fila)\n")
    r_dual = tnep.resolver_bilevel(S, esc, metodo="dualidad")
    r_orac = tnep.resolver_bilevel(S, esc, metodo="oraculo")
    pub = G.PAPER_TABLA_IV["Bilevel"]
    p("| Formulacion | Estado | Plan | Inversion total (M EUR) | Bienestar medio (M EUR/ano) | Objetivo (1) (M EUR/ano) | CPU (s) | max abs(phi) |")
    p("|---|---|---|---:|---:|---:|---:|---:|")
    p(f"| **Paper (Tabla IV)** | - | {pub['plan']} | {pub['inv']:.2f} | {pub['swm']:.2f} | - | {pub['cpu']:.2f} (CPLEX 11) | - |")
    for nombre, r in [("Reproducido: dualidad (36)-(45)", r_dual), ("Reproducido: oraculo (un nivel)", r_orac)]:
        p(f"| {nombre} | {r.estado} | {r.plan_txt()} | {r.inv_total:.3f} | {r.sw_prom:.3f} | {r.obj:.4f} | {r.t:.2f} (Gurobi) | "
          f"{'-' if r.max_phi != r.max_phi else f'{r.max_phi:.1f}'} |")
    assert r_dual.estado == "OPTIMO" and r_orac.estado == "OPTIMO"
    assert abs(r_dual.obj - r_orac.obj) < 1e-6, "dualidad y oraculo deben coincidir"
    assert mismo_plan(r_dual.plan, parse_plan(pub["plan"]))
    assert abs(r_dual.sw_prom - pub["swm"]) < 0.005
    p(f"\nVerificacion: el MILP de dualidad fuerte y el oraculo de un nivel dan el mismo objetivo "
      f"({r_dual.obj:.6f} M EUR/ano). Plan, inversion (25.10) y bienestar medio (292.71) coinciden con el paper.\n")
    out["bilevel"] = dict(plan=r_dual.plan_txt(), inv=r_dual.inv_total, swm=r_dual.sw_prom, obj=r_dual.obj,
                          t=r_dual.t, t_oraculo=r_orac.t, max_phi=r_dual.max_phi)

    # sensibilidad a M, Gamma
    p("### Sensibilidad a las constantes M y Gamma_max (paper: 1000)\n")
    p("| M | Gamma_max | Plan | Objetivo | CPU (s) |")
    p("|---:|---:|---|---:|---:|")
    for Mv, Gv in [(1000, 1000), (10, 1000), (1000, 100000), (10, 100000)]:
        r = tnep.resolver_bilevel(S, esc, metodo="dualidad", M=Mv, gamma=Gv)
        p(f"| {Mv} | {Gv} | {r.plan_txt()} | {r.obj:.6f} | {r.t:.2f} |")
        assert abs(r.obj - r_dual.obj) < 1e-6
    p("\nEl optimo no depende de M ni de Gamma_max en este rango (en el optimo no hay racionamiento y max abs(phi) es pequeno).\n")

    # ------------------------------------------------------------------
    p("## 2. Minimizacion de costos por escenario (Tabla IV, filas 2-4)\n")
    p("Objetivo usado: `0.10*c*x + 8760*cU*r` (inversion anualizada + racionamiento anual), sin presupuesto (2); "
      "es la ponderacion que reproduce las inversiones publicadas (la version sin anualizar no construye lineas).\n")
    filas = []
    for k, sc in enumerate(esc):
        nombre = f"Min-costos Esc. {sc.nombre}"
        pubf = G.PAPER_TABLA_IV[nombre]
        r = tnep.resolver_costmin(S, sc, ponderacion="anual", pool_n=500)
        plan_pub = parse_plan(pubf["plan"])
        en_pool = any(mismo_plan(pl, plan_pub) for pl in r.pool)
        # bienestar medio (3 escenarios) con el plan FIJO: LP primal (no depende de Gamma)
        swm_pub, sws_pub, rac_pub = tnep.evaluar_plan(S, plan_pub, esc)
        swm_mio, sws_mio, rac_mio = tnep.evaluar_plan(S, r.plan, esc)
        sw_claro = S.sigma * sws_pub[k] / 1e6
        filas.append((nombre, pubf, r, en_pool, swm_pub, swm_mio, sw_claro, len(r.pool)))
        out[nombre] = dict(plan_solver=r.plan_txt(), inv=r.inv_total, optimos_empatados=len(r.pool),
                           plan_paper_en_pool=en_pool, swm_plan_paper=swm_pub, swm_plan_solver=swm_mio,
                           sw_mercado_despejado_plan_paper=sw_claro)
        assert abs(r.inv_total - pubf["inv"]) < 0.006, "inversion de costo-min no coincide"
    p("| Aproximacion | Plan paper | Inv. paper | Inv. reprod. | Optimos empatados | Plan paper entre ellos | Bienestar medio paper | Bienestar medio (plan paper, LP primal) | Bienestar medio (plan del solver) |")
    p("|---|---|---:|---:|---:|:--:|---:|---:|---:|")
    for nombre, pubf, r, en_pool, swm_pub, swm_mio, sw_claro, npool in filas:
        # el paper publica 2 decimales: se acepta |dif| <= 0.01 (p. ej. -2613.0261 vs -2613.02)
        marca = "OK" if abs(swm_pub - pubf["swm"]) <= 0.01 else "DIF"
        assert marca == "OK" and en_pool, (nombre, swm_pub, pubf["swm"], en_pool)
        p(f"| {nombre} | {pubf['plan']} | {pubf['inv']:.2f} | {r.inv_total:.3f} | {npool} | {'si' if en_pool else 'NO'} | "
          f"{pubf['swm']:.2f} | {swm_pub:.4f} ({marca}) | {swm_mio:.2f} |")
    p("\nNotas:")
    p("* La inversion minima coincide en los tres escenarios. Hay **optimos alternativos de igual costo** "
      "(2-6 y 4-6 cuestan lo mismo, 5 792.40; 2-3 y 3-5 cuestan lo mismo, 3 861.60). El plan del paper pertenece "
      "al conjunto de optimos empatados y, evaluado con el mercado, reproduce el bienestar medio publicado "
      "(-2 613.02 y -423.57 exactos). El plan que devuelve el solver puede ser otro optimo empatado, con distinto bienestar medio.")
    p("* El bienestar medio se calcula, como en el paper, resolviendo el mercado de cada escenario con el plan fijo. "
      "Se usa el LP primal (5)-(19): con racionamiento los precios sombra superan Gamma_max=1000 y la dualidad fuerte (35) "
      "del MILP seria infactible (ver `verificaciones_modelo.py`).\n")

    p("### Columna 'Social welfare' por escenario de la Tabla IV (275.60 / 291.50 / 286.62)\n")
    p("| Plan (paper) | SW paper (M EUR) | Mercado despejado con el plan fijo (M EUR) | Paper <= despejado |")
    p("|---|---:|---:|:--:|")
    for nombre, pubf, r, en_pool, swm_pub, swm_mio, sw_claro, npool in filas:
        p(f"| {nombre}: {pubf['plan']} | {pubf['sw']:.2f} | {sw_claro:.2f} | {'si' if pubf['sw'] <= sw_claro + 1e-6 else 'NO'} |")
    p("\n**No reproducible por construccion**: esas tres cifras salen del despacho del modelo de minimizacion de costos, "
      "cuyo objetivo no contiene precios; el despacho es degenerado (depende del vertice que devuelva el solver). "
      "Se probaron el mercado despejado, el consumo fijo en el minimo (90 %) y el 100 %, y ninguno da esas cifras. "
      "Todas quedan por debajo del bienestar maximo (mercado despejado), coherente con un despacho no optimo para el mercado. "
      "Las cifras que determinan las conclusiones (plan, inversion, bienestar medio) si se reproducen.\n")

    # ------------------------------------------------------------------
    p("## 3. Modelo basico y centralizado (tesis, Cuadro 6.4)\n")
    rc = tnep.resolver_centralizado(S, esc[0], frac_demanda=1.0)
    swm_c, _, rac_c = tnep.evaluar_plan(S, rc.plan, esc)
    pubc = G.PAPER_TABLA_IV["Basica y centralizada (tesis)"]
    p("| | Plan | Inversion (M EUR) | Bienestar medio (M EUR/ano) |")
    p("|---|---|---:|---:|")
    p(f"| Tesis | {pubc['plan']} | {pubc['inv']:.2f} | {pubc['swm']:.2f} |")
    p(f"| Reproducido (demanda fija = 100 % del escenario base, sin racionamiento) | {rc.plan_txt()} | {rc.inv_total:.3f} | {swm_c:.2f} |")
    assert mismo_plan(rc.plan, parse_plan(pubc["plan"])) and abs(swm_c - pubc["swm"]) < 0.006
    p("\nLa tesis solo dice 'demanda constante igual al pico previsto'; la definicion que reproduce la cifra es el escenario 1 al 100 %.\n")
    out["centralizado"] = dict(plan=rc.plan_txt(), inv=rc.inv_total, swm=swm_c)

    # ------------------------------------------------------------------
    p("## 4. Contingencias simples (Tabla V)\n")
    p("Cuatro escenarios (pesos 0.40/0.25/0.25/0.10); el cuarto saca de servicio el circuito existente indicado.\n")
    p("Se resuelve con la formulacion de dualidad (36)-(45) y, como control, con el oraculo de un nivel (solo primal). "
      "'Racion. max.' es el mayor racionamiento (MW) de un escenario en el plan optimo; 'max abs(phi)' es el mayor |phi| del dual.\n")
    p("| Linea fuera | Plan paper | Plan reproducido | Inv. paper | Inv. reprod. | Bienestar medio paper | Bienestar medio reprod. | Objetivo del plan paper = optimo | Obj. dualidad | Obj. oraculo | Racion. max. (MW) | max abs(phi) | CPU paper (s) | CPU Gurobi dualidad / oraculo (s) |")
    p("|---|---|---|---:|---:|---:|---:|:--:|---:|---:|---:|---:|---:|---:|")
    out["contingencias"] = {}
    for cor, pubf in G.PAPER_TABLA_V.items():
        escc = G.escenarios_contingencia(cor)
        r = tnep.resolver_bilevel(S, escc, metodo="dualidad")
        ro = tnep.resolver_bilevel(S, escc, metodo="oraculo")
        plan_pub = parse_plan(pubf["plan"])
        obj_pub = tnep.objetivo_plan(S, plan_pub, escc)
        equiv = abs(obj_pub - r.obj) < 1e-6
        igual = mismo_plan(r.plan, plan_pub)
        p(f"| {cor[0]}-{cor[1]} | {pubf['plan']} | {r.plan_txt()}{'' if igual else ' (*)'} | {pubf['inv']:.2f} | {r.inv_total:.2f} | "
          f"{pubf['swm']:.2f} | {r.sw_prom:.2f} | {'si' if equiv else 'NO'} | {r.obj:.6f} | {ro.obj:.6f} | "
          f"{max(r.rac_esc):.2f} | {r.max_phi:.1f} | {pubf['cpu']:.2f} | {r.t:.2f} / {ro.t:.2f} |")
        out["contingencias"][f"{cor[0]}-{cor[1]}"] = dict(plan=r.plan_txt(), inv=r.inv_total, swm=r.sw_prom, obj=r.obj,
                                                         obj_oraculo=ro.obj, plan_paper_optimo=equiv, t=r.t,
                                                         racion_max=max(r.rac_esc), max_phi=r.max_phi)
        assert r.estado == "OPTIMO" and ro.estado == "OPTIMO"
        assert ro.obj >= r.obj - 1e-6          # el oraculo relaja el dual: nunca puede ser peor
        if ro.obj > r.obj + 1e-6:
            p(f"\n**HALLAZGO** (linea {cor[0]}-{cor[1]}): el oraculo supera a la formulacion de dualidad "
              f"({ro.obj:.6f} > {r.obj:.6f}); plan del oraculo: {ro.plan_txt()}.\n")
    p("\n(*) plan distinto al publicado: se indica en la columna 'Objetivo del plan paper = optimo' si el plan del paper "
      "alcanza el mismo objetivo (optimo alternativo).\n")

    # ------------------------------------------------------------------
    p("## 5. Objetivo alternativo (20) (sin resultado publicado para comparar)\n")
    r20 = tnep.resolver_bilevel(S, esc, metodo="dualidad", objetivo="racionamiento")
    p(f"Max (20) = -(racionamiento medio) - inversion: plan = {r20.plan_txt()}, inversion = {r20.inv_total:.3f} M EUR, "
      f"bienestar medio = {r20.sw_prom:.2f} M EUR/ano, estado = {r20.estado}.\n")
    out["objetivo_20"] = dict(plan=r20.plan_txt(), inv=r20.inv_total, swm=r20.sw_prom)

    with open(os.path.join(DIR_RES, "garver_resultados.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(MD))
    with open(os.path.join(DIR_RES, "garver_resultados.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    print("\nArchivos escritos en", os.path.abspath(DIR_RES))


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"Tiempo total: {time.time() - t0:.1f} s")
