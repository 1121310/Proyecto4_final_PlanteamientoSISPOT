# -*- coding: utf-8 -*-
"""
RETO 1 - Validacion de los resultados de Garces et al. (2009) en el IEEE RTS de
24 barras: Tabla VIII (sin contingencias), Cuadro 6.9 de la tesis (contingencia
simple) y Tabla X (multicontingencia, 22 escenarios).

Uso:
    python reto1_rts24.py sin                  # Tabla VIII
    python reto1_rts24.py simple               # Cuadro 6.9 (6 corridas)
    python reto1_rts24.py multi [limite_s]     # Tabla X (22 escenarios)
Salida: resultados/rts24_<caso>.md y .json
"""
import json
import os
import sys
import time

import datos_rts24 as R
import tnep
from sistema import parse_plan

DIR_RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resultados")
os.makedirs(DIR_RES, exist_ok=True)
MD = []


def p(texto=""):
    print(texto, flush=True)
    MD.append(texto)


def mismo_plan(a, b):
    return sorted(a.items()) == sorted(b.items())


def guardar(nombre, out):
    with open(os.path.join(DIR_RES, f"rts24_{nombre}.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(MD))
    with open(os.path.join(DIR_RES, f"rts24_{nombre}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)


def fmt_gap(r):
    return "-" if r.gap != r.gap else f"{100 * r.gap:.4f} %"


# ----------------------------------------------------------------------
def caso_sin():
    S = R.sistema_rts24(cmax=23.0, bloques="exactos")
    Si = R.sistema_rts24(cmax=23.0, bloques="impresos")
    esc = R.escenarios_sin_contingencia()
    out = {}
    p("# Reto 1 - IEEE RTS-24, caso sin contingencias (Tabla VIII)\n")
    p(f"Circuitos candidatos (binarios): {len(S.candidatas)}; escenarios: {len(esc)}; presupuesto = {S.cmax} M EUR; "
      f"M = Gamma_max = {S.M:.0f}.\n")
    p("**Dos variantes de los bloques de la Tabla VI** (mismos totales por barra y mismos precios):\n")
    p("* *impresos*: los tamanos de bloque tal como aparecen en la tabla (redondeados, p. ej. 259 / 32.5 / 32.5 en la barra 1).")
    p("* *exactos*: reparto 80/10/10 % de la demanda total y 40/30/30 % de la capacidad de cada generador "
      "(barra 1: 259.2 / 32.4 / 32.4 y 230.4 / 172.8 / 172.8). **Es una inferencia, no un dato publicado**: la tabla impresa es el "
      "redondeo de este reparto (desviacion maxima 0.6 MW) y, como se ve en la seccion 3, reproduce las cuatro cifras de "
      "bienestar medio de la Tabla VIII con error <= 0.003 M EUR, cosa que los bloques impresos no logran.\n")

    pub = R.PAPER_TABLA_VIII["Bilevel"]
    r_d = tnep.resolver_bilevel(S, esc, metodo="dualidad")
    r_o = tnep.resolver_bilevel(S, esc, metodo="oraculo")
    r_i = tnep.resolver_bilevel(Si, esc, metodo="dualidad")
    p("## 1. Modelo bilevel\n")
    p("| Formulacion | Estado | Plan | Inversion (M EUR) | Bienestar medio (M EUR/ano) | Objetivo (1) | CPU (s) | max abs(phi) |")
    p("|---|---|---|---:|---:|---:|---:|---:|")
    p(f"| **Paper (Tabla VIII)** | - | {pub['plan']} | {pub['inv']:.2f} | {pub['swm']:.2f} | - | {pub['cpu']:.2f} (CPLEX 11) | - |")
    for nombre, r in [("Reproducido, bloques exactos: dualidad (36)-(45)", r_d),
                      ("Reproducido, bloques exactos: oraculo", r_o),
                      ("Reproducido, bloques impresos: dualidad", r_i)]:
        p(f"| {nombre} | {r.estado} | {r.plan_txt()} | {r.inv_total:.3f} | {r.sw_prom:.3f} | {r.obj:.4f} | {r.t:.1f} (Gurobi) | "
          f"{'-' if r.max_phi != r.max_phi else f'{r.max_phi:.1f}'} |")
    assert r_d.estado == "OPTIMO" and abs(r_d.obj - r_o.obj) < 1e-6
    assert mismo_plan(r_d.plan, parse_plan(pub["plan"])) and mismo_plan(r_i.plan, parse_plan(pub["plan"]))
    assert abs(r_d.sw_prom - pub["swm"]) < 0.005
    out["bilevel"] = dict(plan=r_d.plan_txt(), inv=r_d.inv_total, swm=r_d.sw_prom, obj=r_d.obj, t=r_d.t,
                          t_oraculo=r_o.t, max_phi=r_d.max_phi,
                          impresos=dict(plan=r_i.plan_txt(), swm=r_i.sw_prom, obj=r_i.obj, t=r_i.t))
    p("\nCon ambas variantes se obtiene el plan publicado (6-10, 10-11, 11-13, 14-23) y su inversion (22.99). "
      "El bienestar medio con bloques exactos coincide con el paper (742.59).\n")

    p("## 2. Minimizacion de costos por escenario\n")
    p("Mismo objetivo que en Garver: `0.10*c*x + 8760*cU*r`, sin presupuesto.\n")
    p("| Aproximacion | Plan paper | Inv. paper | Plan del solver | Inv. reprod. | Optimos empatados | Plan paper entre ellos | CPU (s) |")
    p("|---|---|---:|---|---:|---:|:--:|---:|")
    out["costmin"] = {}
    for k, sc in enumerate(esc):
        nombre = f"Min-costos Esc. {sc.nombre}"
        pf = R.PAPER_TABLA_VIII[nombre]
        r = tnep.resolver_costmin(S, sc, ponderacion="anual", pool_n=200, time_limit=600)
        en_pool = any(mismo_plan(pl, parse_plan(pf["plan"])) for pl in r.pool)
        p(f"| {nombre} | {pf['plan']} | {pf['inv']:.2f} | {r.plan_txt()} | {r.inv_total:.3f} | {len(r.pool)} | "
          f"{'si' if en_pool else 'NO'} | {r.t:.1f} |")
        out["costmin"][nombre] = dict(plan_solver=r.plan_txt(), inv=r.inv_total, optimos_empatados=len(r.pool),
                                      plan_paper_en_pool=en_pool, estado=r.estado)
        assert abs(r.inv_total - pf["inv"]) < 0.006, (nombre, r.inv_total, pf["inv"])

    p("\n## 3. Bienestar medio de cada plan publicado (mercado de los 3 escenarios con el plan fijo)\n")
    p("| Plan | Paper | Bloques impresos | Dif. | Bloques exactos | Dif. | Racionamiento por escenario (MW) |")
    p("|---|---:|---:|---:|---:|---:|---|")
    dif_i, dif_e = [], []
    out["swm_planes"] = {}
    for nombre, pf in R.PAPER_TABLA_VIII.items():
        swi, _, _ = tnep.evaluar_plan(Si, parse_plan(pf["plan"]), esc)
        swe, _, racs = tnep.evaluar_plan(S, parse_plan(pf["plan"]), esc)
        dif_i.append(swi - pf["swm"])
        dif_e.append(swe - pf["swm"])
        p(f"| {nombre}: {pf['plan']} | {pf['swm']:.2f} | {swi:.3f} | {swi - pf['swm']:+.3f} | {swe:.3f} | "
          f"{swe - pf['swm']:+.3f} | {', '.join(f'{x:.1f}' for x in racs)} |")
        out["swm_planes"][nombre] = dict(paper=pf["swm"], impresos=swi, exactos=swe, racion=racs)
    assert max(abs(d) for d in dif_e) < 0.0035, dif_e
    p(f"\nCon los bloques impresos las cuatro cifras quedan desplazadas entre {min(dif_i):+.3f} y {max(dif_i):+.3f} M EUR/ano "
      f"(~0.025 %), un desplazamiento casi constante, igual en planes con y sin racionamiento (no es un efecto de red). "
      f"Con los bloques exactos el error maximo es {max(abs(d) for d in dif_e):.3f} M EUR/ano (<= redondeo del paper). "
      "Conclusion: el modelo reproduce el paper; la discrepancia con los datos impresos es de redondeo de la Tabla VI. "
      "Ningun parametro del modelo se ajusto para forzar el acuerdo.\n")
    guardar("sin_contingencia", out)


# ----------------------------------------------------------------------
def caso_simple():
    S = R.sistema_rts24(cmax=23.0, bloques="exactos")
    Si = R.sistema_rts24(cmax=23.0, bloques="impresos")
    out = {}
    p("# Reto 1 - IEEE RTS-24, contingencia simple (tesis, Cuadro 6.9)\n")
    p("Cuatro escenarios (pesos 0.40/0.25/0.25/0.10); presupuesto 23 M EUR. El paper no publica este caso; los datos de "
      "comparacion son del Cuadro 6.9 de la tesis. Se corre con bloques *exactos* (inferencia, ver reporte del caso sin "
      "contingencias) y con bloques *impresos* de la Tabla VI; ademas de la formulacion de dualidad se corre el oraculo de "
      "un nivel como control.\n")
    p("| Linea fuera | Bloques | Plan tesis | Plan reproducido | Inv. tesis | Inv. reprod. | Bienestar medio tesis | Bienestar medio reprod. | Objetivo del plan tesis = optimo | Obj. dualidad | Obj. oraculo | Racion. max. (MW) | max abs(phi) | CPU tesis (min) | CPU Gurobi dualidad / oraculo (min) |")
    p("|---|---|---|---|---:|---:|---:|---:|:--:|---:|---:|---:|---:|---:|---:|")
    for cor, pf in R.PAPER_CUADRO_6_9.items():
        escc = R.escenarios_contingencia_simple(cor)
        for etiqueta, Sx in (("exactos", S), ("impresos", Si)):
            r = tnep.resolver_bilevel(Sx, escc, metodo="dualidad", time_limit=3600)
            ro = tnep.resolver_bilevel(Sx, escc, metodo="oraculo", time_limit=3600)
            plan_pub = parse_plan(pf["plan"])
            obj_pub = tnep.objetivo_plan(Sx, plan_pub, escc)
            equiv = abs(obj_pub - r.obj) < 1e-6
            igual = mismo_plan(r.plan, plan_pub)
            p(f"| {cor[0]}-{cor[1]} | {etiqueta} | {pf['plan']} | {r.plan_txt()}{'' if igual else ' (*)'} | {pf['inv']:.2f} | {r.inv_total:.2f} | "
              f"{pf['swm']:.2f} | {r.sw_prom:.2f} | {'si' if equiv else 'NO'} | {r.obj:.6f} | {ro.obj:.6f} | "
              f"{max(r.rac_esc):.2f} | {r.max_phi:.1f} | {pf['cpu_min']:.2f} | {r.t / 60:.2f} / {ro.t / 60:.2f} |")
            out[f"{cor[0]}-{cor[1]}_{etiqueta}"] = dict(
                estado=r.estado, plan=r.plan_txt(), inv=r.inv_total, swm=r.sw_prom, obj=r.obj, obj_oraculo=ro.obj,
                obj_plan_tesis=obj_pub, plan_tesis_optimo=equiv, t=r.t, t_oraculo=ro.t, gap=r.gap,
                racion_max=max(r.rac_esc), max_phi=r.max_phi)
            assert ro.obj >= r.obj - 1e-6
            if ro.obj > r.obj + 1e-6:
                p(f"\n**HALLAZGO** (linea {cor[0]}-{cor[1]}, bloques {etiqueta}): el oraculo supera a la formulacion de "
                  f"dualidad ({ro.obj:.6f} > {r.obj:.6f}); plan del oraculo: {ro.plan_txt()}.\n")
    p("\n(*) plan distinto al de la tesis (solo en el orden de escritura o por optimo alternativo); la columna "
      "'Objetivo del plan tesis = optimo' indica si el plan de la tesis alcanza el mismo objetivo.\n")
    guardar("contingencia_simple", out)


# ----------------------------------------------------------------------
def caso_multi(limite=3600):
    from sistema import Escenario
    S = R.sistema_rts24(cmax=30.0)
    esc_imp = R.escenarios_multicontingencia()
    # pesos 'exactos': 0.25 y 0.75/21 = 1/28 (suma 1). Los '0.0357' impresos son 1/28 redondeado.
    esc_ex = [Escenario(e.nombre, 0.25 if i == 0 else 0.75 / 21, e.factores, e.linea_fuera)
              for i, e in enumerate(esc_imp)]
    pub = R.PAPER_TABLA_X
    out = {}
    p("# Reto 1 - IEEE RTS-24, multicontingencia (Tabla X)\n")
    p(f"22 escenarios (16 de demanda + 6 contingencias), presupuesto 30 M EUR. Limite de tiempo por corrida: {limite} s.\n")
    p(f"**Dos variantes de pesos**: *impresos* (0.25 y 21 x 0.0357, suma = {sum(e.peso for e in esc_imp):.4f}) y *exactos* "
      f"(0.25 y 0.75/21 = 1/28 = 0.035714, suma = 1). La segunda es una inferencia: los 0.0357 impresos son 1/28 redondeado; "
      "se valida mas abajo porque reproduce el bienestar medio de la Tabla X (717.97) y la suma de probabilidades debe ser 1.\n")
    plan_pub = parse_plan(pub["plan"])
    sw_imp = tnep.evaluar_plan(S, plan_pub, esc_imp)[0]
    sw_ex = tnep.evaluar_plan(S, plan_pub, esc_ex)[0]
    sw_imp_bi = tnep.evaluar_plan(R.sistema_rts24(cmax=30.0, bloques="impresos"), plan_pub, esc_imp)[0]
    sw_ex_bi = tnep.evaluar_plan(R.sistema_rts24(cmax=30.0, bloques="impresos"), plan_pub, esc_ex)[0]
    p(f"Bienestar medio del plan publicado evaluado con este modelo (paper = {pub['swm']:.2f}):\n")
    p("| Bloques | Pesos impresos | Pesos exactos |")
    p("|---|---:|---:|")
    p(f"| exactos | {sw_imp:.3f} | {sw_ex:.3f} |")
    p(f"| impresos | {sw_imp_bi:.3f} | {sw_ex_bi:.3f} |")
    p("")
    out["swm_plan_paper"] = dict(paper=pub["swm"], bloques_exactos=dict(pesos_impresos=sw_imp, pesos_exactos=sw_ex),
                                 bloques_impresos=dict(pesos_impresos=sw_imp_bi, pesos_exactos=sw_ex_bi))
    Si = R.sistema_rts24(cmax=30.0, bloques="impresos")
    p("Bloques: 'exactos' = inferencia (80/10/10 % y 40/30/30 %); 'impresos' = Tabla VI tal cual. Se reportan ambas.\n")
    p("| Formulacion | Bloques | Pesos | Estado | Plan | Inversion (M EUR) | Bienestar medio (M EUR/ano) | Objetivo | Gap | CPU |")
    p("|---|---|---|---|---|---:|---:|---:|---:|---:|")
    p(f"| **Paper (Tabla X)** | - | - | - | {pub['plan']} | {pub['inv']:.2f} | {pub['swm']:.2f} | - | - | {pub['cpu_h']:.2f} h (CPLEX 11) |")
    corridas = [("oraculo (un nivel)", "oraculo", "impresos", "impresos", Si, esc_imp),
                ("oraculo (un nivel)", "oraculo", "exactos", "impresos", S, esc_imp),
                ("oraculo (un nivel)", "oraculo", "exactos", "exactos", S, esc_ex),
                ("dualidad (36)-(45)", "dualidad", "exactos", "exactos", S, esc_ex)]
    for nombre, metodo, bl, pesos, Sx, esc in corridas:
        r = tnep.resolver_bilevel(Sx, esc, metodo=metodo, time_limit=limite)
        p(f"| Reproducido: {nombre} | {bl} | {pesos} | {r.estado} | {r.plan_txt()} | {r.inv_total:.3f} | {r.sw_prom:.3f} | {r.obj:.4f} | "
          f"{fmt_gap(r)} | {r.t:.1f} s |")
        out[f"{metodo}_bloques-{bl}_pesos-{pesos}"] = dict(
            estado=r.estado, plan=r.plan_txt(), inv=r.inv_total, swm=r.sw_prom, obj=r.obj, gap=r.gap, t=r.t,
            max_phi=r.max_phi, racion_max=max(r.rac_esc) if r.rac_esc else None)
        guardar("multicontingencia", out)
    p("")
    guardar("multicontingencia", out)


if __name__ == "__main__":
    caso = sys.argv[1] if len(sys.argv) > 1 else "sin"
    t0 = time.time()
    if caso == "sin":
        caso_sin()
    elif caso == "simple":
        caso_simple()
    elif caso == "multi":
        caso_multi(float(sys.argv[2]) if len(sys.argv) > 2 else 3600)
    else:
        raise SystemExit("uso: reto1_rts24.py sin|simple|multi [limite_s]")
    print(f"\nTiempo total: {time.time() - t0:.1f} s")
