# -*- coding: utf-8 -*-
"""
Reto 2 - incertidumbre en el TNEP bilevel; validacion de Heydari y Barforoushi
(EPSR 229, 2024) para Garver (Tabla 4, Figs. 4-5).

  python reto2_garver.py seg     # metodo 1: curva de duracion + 6 vientos/bloque (el del paper)
  python reto2_garver.py dual    # los 4 casos con la formulacion de dualidad (55)-(70)
  python reto2_garver.py pem     # metodo 2: Point Estimate Method de Hong (1998)
  python reto2_garver.py mc      # metodo 3: Monte Carlo (Beta por bloque) + reduccion k-means
  python reto2_garver.py sens    # sensibilidad al coeficiente de demanda (Figs. 4-5)
  python reto2_garver.py all
Salidas: resultados/reto2_*.md y .json
"""
import json
import sys
import time
from pathlib import Path

import datos_heydari as dh
from tnep_heydari import (resolver_h, evaluar_h, rango_energia_h, parse_plan, plan_txt, costo_plan)
from incertidumbre import segmentacion, pem_casos, montecarlo_kmeans

RES = Path(__file__).resolve().parent.parent / "resultados"
RES.mkdir(exist_ok=True)
E_EOL_PAPER = 2213.81            # GWh, Tabla 4 (casos 2 y 3)


def _viento(escala=1.0):
    return [[w * escala for w in b] for b in dh.VIENTO]


def preparar(caso, escala_viento=1.0, **kw):
    S = dh.sistema_heydari(caso, **kw)
    casos, grupos = segmentacion(dh.N_H, dh.CARGA, _viento(escala_viento), con_viento=(caso != 1))
    beta = dh.BETA if caso == 4 else 0.0
    return S, casos, grupos, beta


def pct(a, b):
    return f"{100.0 * (a / b - 1.0):+.2f} %" if b else "-"


def guardar(nombre, md, datos):
    (RES / f"{nombre}.md").write_text(md, encoding="utf-8")
    (RES / f"{nombre}.json").write_text(json.dumps(datos, indent=1, default=str), encoding="utf-8")
    print(f"\n-> resultados/{nombre}.md")


# ----------------------------------------------------------------------
# 1) segmentacion por curva de duracion (metodo del paper)
# ----------------------------------------------------------------------
def seccion_seg():
    out, md = {}, []
    md.append("# Reto 2 - Garver de Heydari (2024): metodo 1, segmentacion por curva de duracion\n")
    md.append("5 bloques de la curva de duracion (1000/3000/3000/1000/760 h) x 6 intensidades de viento "
              "equiprobables = 30 mercados (caso 1, sin viento: 5). Carga y viento **digitalizados de la "
              "Fig. 2** (datos inferidos, ~0.002 p.u.). Comparacion con la Tabla 4 en dos pruebas por caso: "
              "(a) *prueba de datos*: el plan publicado evaluado con mis datos; (b) *prueba de optimalidad*: "
              "el optimo de mi modelo (MILP de un nivel, verificado contra la formulacion de dualidad en "
              "`reto2_garver.py dual`) frente al plan publicado.\n")
    e_dig = sum(nh * sum(w) / len(w) for nh, w in zip(dh.N_H, dh.VIENTO)) * 1050.0 / 1e3
    k_viento = E_EOL_PAPER / e_dig
    for caso in (1, 2, 3, 4):
        ref = dh.PAPER_TABLA_4[caso]
        S, casos, grupos, beta = preparar(caso)
        pa, pd = parse_plan(ref["ac"]), parse_plan(ref["dc"])
        ev = evaluar_h(S, pa, pd, casos, grupos, dh.ALFA, beta, detalle=True)
        rg = rango_energia_h(S, pa, pd, casos)
        opt = resolver_h(S, casos, grupos, dh.ALFA, beta, metodo="oraculo")
        evo = evaluar_h(S, opt.plan_ac, opt.plan_dc, casos, grupos, dh.ALFA, beta)
        mismo = (opt.plan_ac == pa and opt.plan_dc == pd)
        r = dict(caso=caso, paper=ref, ev={k: v for k, v in ev.items() if k != "det"}, rango=rg,
                 optimo=dict(estado=opt.estado, plan=opt.plan_txt(), inv=opt.inv_total, esw=opt.esw,
                             obj=opt.obj, cvar=opt.cvar, t=opt.t, gap=opt.gap),
                 mismo_plan=mismo, regret=opt.obj - ev["obj"], consistencia=opt.obj - evo["obj"])
        # variante: 4 HVAC por corredor contando el existente
        S4, c4, g4, b4 = preparar(caso, ac_total=True)
        o4 = resolver_h(S4, c4, g4, dh.ALFA, b4, metodo="oraculo")
        r["ac_total"] = dict(plan=o4.plan_txt(), obj=o4.obj, igual=(o4.plan_ac == opt.plan_ac and o4.plan_dc == opt.plan_dc))
        # sensibilidad (NO es ajuste): viento escalado por el factor implicito en la energia eolica publicada
        if caso >= 2:
            Sk, ck, gk, bk = preparar(caso, escala_viento=k_viento)
            evk = evaluar_h(Sk, pa, pd, ck, gk, dh.ALFA, bk, detalle=True)
            ok = resolver_h(Sk, ck, gk, dh.ALFA, bk, metodo="oraculo")
            r["viento_escalado"] = dict(k=k_viento, esw=evk["esw"], obj=evk["obj"], cvar=evk["cvar"],
                                        eolica=evk["GWh_eolica"], plan=ok.plan_txt(),
                                        mismo_plan=(ok.plan_ac == pa and ok.plan_dc == pd))
        out[caso] = r

        md.append(f"\n## Caso {caso}" + {1: " (solo termicas)", 2: " (eolica, solo HVAC)",
                                          3: " (eolica, HVAC + HVDC)", 4: " (eolica, HVAC + HVDC, CVaR)"}[caso])
        md.append("\n| | Paper (Tabla 4) | Reproducido, plan publicado | Dif. |\n|---|---:|---:|---:|")
        md.append(f"| Plan HVAC / HVDC | {ref['ac']} / {ref['dc'] or '-'} | {plan_txt(pa)} / {plan_txt(pd)} | |")
        md.append(f"| Inversion total (M$) | {ref['inv']} | {ev['inv_total']:.1f} | exacta |")
        md.append(f"| Bienestar esperado ESW (M$) | {ref['esw']} | {ev['esw']:.3f} | {pct(ev['esw'], ref['esw'])} |")
        if beta:
            md.append(f"| CVaR (M$) | {ref['cvar']} | {ev['cvar']:.3f} | {pct(ev['cvar'], ref['cvar'])} |")
        md.append(f"| Funcion objetivo del nivel superior (M$) | {ref['obj']} | {ev['obj']:.3f} | {pct(ev['obj'], ref['obj'])} |")
        def rango_txt(v, lo, hi):
            if lo - 0.05 <= v <= hi + 0.05:
                return "dentro del rango"
            cerca = lo if v < lo else hi
            return f"fuera del rango ({pct(cerca, v)} el limite mas cercano)"
        md.append(f"| Energia termica (GWh) | {ref['fosil']} | {ev['GWh_termica']:.1f} (rango de empates {rg['termica'][0]:.1f}-{rg['termica'][1]:.1f}) | "
                  + rango_txt(ref['fosil'], *rg['termica']) + " |")
        md.append(f"| Energia eolica (GWh) | {ref['eolica'] or '-'} | {ev['GWh_eolica']:.1f} | {pct(ev['GWh_eolica'], ref['eolica']) if ref['eolica'] else '-'} |")
        md.append(f"| Consumo (GWh) | {ref['consumo']} | {ev['GWh_consumo']:.1f} (rango {rg['consumo'][0]:.1f}-{rg['consumo'][1]:.1f}) | "
                  + rango_txt(ref['consumo'], *rg['consumo']) + " |")
        if caso >= 2:
            pub_rac = ref['consumo'] - ref['fosil'] - ref['eolica']
            md.append(f"| Racionamiento (GWh) = consumo - termica - eolica | {pub_rac:.2f} | {ev['GWh_racion']:.2f} | {pct(ev['GWh_racion'], pub_rac)} |")
        md.append(f"| LMP medio ($/MWh) | {ref['lmp']} | {ev['lmp_medio']:.3f} | {pct(ev['lmp_medio'], ref['lmp'])} |")
        md.append(f"| Desv. tipica de LMP ($/MWh) | {ref['lmp_std']} | {ev['lmp_std']:.3f} | {pct(ev['lmp_std'], ref['lmp_std'])} |")
        md.append(f"\n*Prueba de optimalidad.* Optimo de mi modelo: **{opt.plan_txt()}**, inversion {opt.inv_total:.1f} M$, "
                  f"objetivo {opt.obj:.3f} M$ ({opt.estado}, gap {opt.gap:.1e}, {opt.t:.1f} s; recalculado con LP de plan fijo: "
                  f"diferencia {opt.obj - evo['obj']:.1e}). "
                  + ("**Coincide con el plan publicado** (arrepentimiento 0)." if mismo else
                     f"Difiere del publicado; con mis datos el plan publicado queda {opt.obj - ev['obj']:.3f} M$ por debajo del optimo."))
        md.append(f"\nVariante *4 circuitos HVAC por corredor contando el existente*: plan {o4.plan_txt()} "
                  + ("(igual)." if r["ac_total"]["igual"] else "(DIFERENTE)."))
        if caso >= 2:
            v = r["viento_escalado"]
            md.append(f"\n*Sensibilidad (no es un ajuste)*: con el viento digitalizado escalado por {v['k']:.4f} "
                      f"(el factor que iguala la energia eolica publicada, {E_EOL_PAPER} GWh), el plan publicado da ESW = "
                      f"{v['esw']:.3f} M$ ({pct(v['esw'], ref['esw'])} vs paper)" +
                      (f", CVaR = {v['cvar']:.3f} M$ ({pct(v['cvar'], ref['cvar'])})" if beta else "") +
                      f" y el optimo {'sigue siendo el publicado' if v['mismo_plan'] else 'cambia a ' + v['plan']}.")
        print(f"caso {caso}: mismo plan={mismo}  ESW {ev['esw']:.3f} (paper {ref['esw']})  obj {ev['obj']:.3f} (paper {ref['obj']})")
    md.append("\n## Notas\n")
    md.append("* **LMP**: el paper no define sus estadisticos. La media y la desviacion tipica poblacional de los "
              "pares (caso, barra), sin ponderar por horas, reproducen la Tabla 4 en los casos 1-3 con 4-6 cifras; "
              "en el caso 4 difiere 0.2-0.3 %. Los LMP son constantes a trozos, por eso no dependen de la "
              "resolucion de la digitalizacion; validan precios, red y estructura de racionamiento, no los niveles de viento.")
    md.append("* **Energias**: hay empates de precio (p. ej. 77 $/MWh en una oferta y una demanda), de modo que el "
              "despacho optimo no es unico. Se da el rango sobre todos los despachos igualmente optimos; los valores "
              "publicados de energia termica y consumo caen dentro. El racionamiento y la energia eolica si son unicos.")
    guardar("reto2_segmentacion", "\n".join(md), out)
    return out


# ----------------------------------------------------------------------
# 2) formulacion de dualidad del paper en los 4 casos
# ----------------------------------------------------------------------
def seccion_dual(limite=3600.0):
    out, md = {}, ["# Reto 2 - formulacion de dualidad (55)-(70) frente al MILP de un nivel\n",
                   f"Gamma = Gamma' = {dh.sistema_heydari(2).Gamma:g} (Heydari no da su valor), M = 1000, limite {limite:g} s. "
                   "Se anade una corrida con Gamma = 1000 (el valor que usa Garces) para ver si esa cota es activa: "
                   "con racionamiento los precios sombra llegan a ~10 x el precio de la demanda y los phi a varios miles.\n",
                   "| Caso | Estado (Gamma 1e4) | Plan | Objetivo dualidad | Objetivo un nivel | max abs(phi) | Tiempo dualidad (s) | Con Gamma = 1000 |\n|---|---|---|---:|---:|---:|---:|---|"]
    for caso in (1, 2, 3, 4):
        S, casos, grupos, beta = preparar(caso)
        o = resolver_h(S, casos, grupos, dh.ALFA, beta, metodo="oraculo")
        d = resolver_h(S, casos, grupos, dh.ALFA, beta, metodo="dualidad", time_limit=limite)
        d1 = resolver_h(S, casos, grupos, dh.ALFA, beta, metodo="dualidad", time_limit=limite, gamma=1000.0)
        g1 = (f"{d1.estado}" if d1.estado != "OPTIMO" else
              f"OPTIMO, obj {d1.obj:.4f}" + (" (igual)" if abs(d1.obj - o.obj) < 1e-6 else " (DIFIERE)"))
        # causa: con el plan publicado FIJO, Gamma = 1000 puede hacer infactible la dualidad fuerte (54)
        ref = dh.PAPER_TABLA_4[caso]
        fix = (parse_plan(ref["ac"]), parse_plan(ref["dc"]))
        df = resolver_h(S, casos, grupos, dh.ALFA, beta, metodo="dualidad", fix=fix, time_limit=limite, gamma=1000.0)
        gf = df.estado if df.estado != "OPTIMO" else "factible"
        out[caso] = dict(dual=dict(estado=d.estado, plan=d.plan_txt(), obj=d.obj, t=d.t, gap=d.gap, max_phi=d.max_phi),
                         oraculo=dict(plan=o.plan_txt(), obj=o.obj, t=o.t),
                         gamma1000=dict(estado=d1.estado, obj=d1.obj, plan=d1.plan_txt(), plan_publicado_fijo=df.estado))
        md.append(f"| {caso} | {d.estado} | {d.plan_txt()} | {d.obj:.4f} | {o.obj:.4f} | {d.max_phi:.0f} | {d.t:.1f} | {g1}; plan publicado fijo: {gf} |")
        print(f"caso {caso}: dual {d.estado} obj {d.obj:.4f} t={d.t:.1f}s | oraculo {o.obj:.4f} t={o.t:.1f}s | max phi {d.max_phi:.0f} | gamma1000: {g1}; plan publicado fijo: {gf}")
    guardar("reto2_dualidad", "\n".join(md), out)
    return out


# ----------------------------------------------------------------------
# 3) Point Estimate Method de Hong
# ----------------------------------------------------------------------
def seccion_pem():
    out, md = {}, ["# Reto 2 - metodo 2: Point Estimate Method (Hong, 1998)\n",
                   "Para cada bloque de la curva de duracion la intensidad eolica es una variable aleatoria (m = 1) "
                   "cuya distribucion son los 6 niveles equiprobables del paper; el PEM la sustituye por 2 puntos "
                   "(esquema 2m) o 3 puntos (2m+1) con pesos de Hong (ver `incertidumbre.py`). Solo casos 2 y 3: el caso 1 "
                   "no tiene viento y en el caso 4 los pesos del PEM dependen del bloque, de modo que no definen "
                   "escenarios anuales con probabilidad propia para el CVaR.\n"]
    planes_pub = {c: (parse_plan(dh.PAPER_TABLA_4[c]["ac"]), parse_plan(dh.PAPER_TABLA_4[c]["dc"])) for c in (1, 2, 3, 4)}
    S3 = dh.sistema_heydari(3)
    seg3, _ = segmentacion(dh.N_H, dh.CARGA, dh.VIENTO, True)
    # puntos por bloque
    _, info = pem_casos(dh.N_H, dh.CARGA, dh.VIENTO, "2m+1")
    md.append("## Distribucion de cada bloque y puntos del PEM 2m+1\n")
    md.append("| Bloque | media | desv. | asimetria | curtosis | puntos (viento p.u. : peso) |\n|---|---:|---:|---:|---:|---|")
    for inf in info:
        mu, sd, l3, l4 = inf["momentos"]
        md.append(f"| {inf['bloque']} | {mu:.4f} | {sd:.4f} | {l3:.3f} | {l4:.3f} | " +
                  "; ".join(f"{x:.4f}:{p:.3f}" for x, p in inf["puntos"]) + " |")
    out["info"] = [dict(bloque=i["bloque"], momentos=i["momentos"], puntos=i["puntos"]) for i in info]
    # exactitud del estimador de la esperanza para planes fijos
    md.append("\n## Exactitud del estimador del bienestar esperado (ESW, M$) para planes fijos\n")
    md.append("Se evaluan los 4 planes publicados en el sistema con eolica y HVDC (caso 3). El plan del caso 1 "
              "(sin HVDC ni refuerzos para el viento) es insuficiente en ese sistema y racionan mucho: da un ESW "
              "negativo y sirve de prueba de estres. El error es (PEM - 6 niveles) y su % se calcula sobre |ESW de 6 niveles|.\n")
    md.append("| Plan (publicado en caso) | 6 niveles (30 LP) | PEM 2m (10 LP) | error | PEM 2m+1 (15 LP) | error |\n|---|---:|---:|---:|---:|---:|")
    esq_casos = {e: pem_casos(dh.N_H, dh.CARGA, dh.VIENTO, e)[0] for e in ("2m", "2m+1")}
    est = {}
    for c in (1, 2, 3, 4):
        pa, pd = planes_pub[c]
        e6 = evaluar_h(S3, pa, pd, seg3)["esw"]
        fila = [e6]
        for e in ("2m", "2m+1"):
            ep = evaluar_h(S3, pa, pd, esq_casos[e])["esw"]
            fila += [ep, ep - e6, 100 * (ep - e6) / abs(e6)]
        est[c] = fila
        md.append(f"| caso {c}: {dh.PAPER_TABLA_4[c]['ac']} / {dh.PAPER_TABLA_4[c]['dc'] or '-'} | {fila[0]:.3f} | {fila[1]:.3f} | {fila[2]:+.3f} ({fila[3]:+.2f} %) | {fila[4]:.3f} | {fila[5]:+.3f} ({fila[6]:+.2f} %) |")
    out["estimador"] = est
    # optimizacion con cada metodo
    md.append("\n## Plan optimo con cada metodo de incertidumbre\n")
    md.append("El objetivo de la ultima columna es el del plan evaluado con los 30 mercados de referencia (6 niveles), "
              "para comparar planes en la misma vara.\n")
    md.append("| Caso | Metodo | mercados | Estado | Plan HVAC / HVDC | Inv. (M$) | Objetivo propio | Objetivo con 6 niveles | Arrepentimiento (M$) | t (s) |\n|---|---|---:|---|---|---:|---:|---:|---:|---:|")
    out["opt"] = {}
    for caso in (2, 3):
        S, seg, grupos, beta = preparar(caso)
        base = resolver_h(S, seg, grupos, dh.ALFA, 0.0, metodo="oraculo")
        base_obj = evaluar_h(S, base.plan_ac, base.plan_dc, seg)["obj"]
        md.append(f"| {caso} | 6 niveles (paper) | {len(seg)} | {base.estado} | {plan_txt(base.plan_ac)} / {plan_txt(base.plan_dc)} | {base.inv_total:.1f} | {base.obj:.3f} | {base_obj:.3f} | 0 | {base.t:.1f} |")
        out["opt"][f"{caso}-6"] = dict(plan=base.plan_txt(), obj=base.obj, t=base.t)
        for e in ("2m", "2m+1"):
            cp = esq_casos[e]
            r = resolver_h(S, cp, None, dh.ALFA, 0.0, metodo="oraculo")
            ref_obj = evaluar_h(S, r.plan_ac, r.plan_dc, seg)["obj"]
            md.append(f"| {caso} | PEM {e} | {len(cp)} | {r.estado} | {plan_txt(r.plan_ac)} / {plan_txt(r.plan_dc)} | {r.inv_total:.1f} | {r.obj:.3f} | {ref_obj:.3f} | {base_obj - ref_obj:.3f} | {r.t:.1f} |")
            out["opt"][f"{caso}-{e}"] = dict(plan=r.plan_txt(), obj=r.obj, obj_ref=ref_obj, regret=base_obj - ref_obj, t=r.t)
            print(f"caso {caso} PEM {e}: {r.plan_txt()}  regret {base_obj - ref_obj:.3f}  t={r.t:.1f}s")
    guardar("reto2_pem", "\n".join(md), out)
    return out


# ----------------------------------------------------------------------
# 4) sensibilidad al coeficiente de demanda (Figs. 4 y 5)
# ----------------------------------------------------------------------
REL = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]


def seccion_sens():
    out, md = {}, ["# Reto 2 - sensibilidad al coeficiente de demanda (Figs. 4 y 5 de Heydari)\n",
                   "Coeficiente relativo r = 0.85...1.15 aplicado sobre el coeficiente base 2.2 (r = 1 reproduce la "
                   "Tabla 4). Las figuras del paper son graficos de barras 3D sin valores tabulados: se comparan las "
                   "tendencias descritas en el texto, no cifras.\n"]
    for caso in (1, 2, 3, 4):
        out[caso] = []
        for rel in REL:
            S, casos, grupos, beta = preparar(caso, coef_dem=2.2 * rel)
            o = resolver_h(S, casos, grupos, dh.ALFA, beta, metodo="oraculo")
            ev = evaluar_h(S, o.plan_ac, o.plan_dc, casos, grupos, dh.ALFA, beta, detalle=True)
            out[caso].append(dict(rel=rel, plan=o.plan_txt(), inv=o.inv_total, esw=ev["esw"], racion=ev["GWh_racion"], t=o.t))
            print(f"caso {caso} r={rel:.2f}: inv {o.inv_total:.1f} esw {ev['esw']:.2f} rac {ev['GWh_racion']:.1f} | {o.plan_txt()}")
    md.append("| r | " + " | ".join(f"Caso {c}: inv. / ESW" for c in (1, 2, 3, 4)) + " |\n|---|" + "---:|" * 4)
    for i, rel in enumerate(REL):
        md.append(f"| {rel:.2f} | " + " | ".join(f"{out[c][i]['inv']:.1f} / {out[c][i]['esw']:.1f}" for c in (1, 2, 3, 4)) + " |")
    md.append("\n## Planes\n\n| r | " + " | ".join(f"Caso {c}" for c in (1, 2, 3, 4)) + " |\n|---|---|---|---|---|")
    for i, rel in enumerate(REL):
        md.append(f"| {rel:.2f} | " + " | ".join(out[c][i]["plan"] for c in (1, 2, 3, 4)) + " |")
    inv = {c: [x["inv"] for x in out[c]] for c in out}
    esw = {c: [x["esw"] for x in out[c]] for c in out}
    mono = lambda v: all(b >= a - 1e-9 for a, b in zip(v, v[1:]))
    # primer indice desde el cual la inversion ya no cambia (>= 2 puntos), o None
    meseta = lambda v: next((i for i in range(len(v) - 1) if max(v[i:]) - min(v[i:]) < 1e-9), None)
    sube_baja = lambda v: max(v) > v[0] and v[-1] < max(v) - 1e-9
    chk = [
        ("Caso 1: la inversion crece con el coeficiente de demanda", mono(inv[1])),
        ("Caso 1: el bienestar crece uniformemente con el coeficiente", mono(esw[1])),
        ("Casos 2-4: la inversion crece con el coeficiente (hasta el limite de presupuesto)",
         all(mono(inv[c]) for c in (2, 3, 4))),
        ("Casos 2-4: el bienestar primero crece y luego decrece (mas racionamiento)",
         all(sube_baja(esw[c]) for c in (2, 3, 4))),
        ("Casos 2 y 3: misma inversion donde no se instala HVDC (coef. bajos)",
         any(abs(a - b) < 1e-9 for a, b in zip(inv[2], inv[3]))),
        ("La inversion del caso 3 llega a su meseta (limite de presupuesto) con un coeficiente mayor o igual "
         "que la del caso 2, o no llega",
         (lambda f2, f3: f3 is None or (f2 is not None and f3 >= f2))(meseta(inv[2]), meseta(inv[3]))),
        ("Caso 4: invierte menos o igual que el caso 3 en todos los coeficientes",
         all(a <= b + 1e-9 for a, b in zip(inv[4], inv[3]))),
    ]
    md.append("\n## Tendencias descritas en el paper (Sec. 4.1.3)\n")
    for txt, v in chk:
        md.append(f"* {'CONSISTENTE' if v else 'NO CONSISTENTE'}: {txt}")
        print(("OK  " if v else "NO  ") + txt)
    out["chequeos"] = [dict(texto=t, ok=bool(v)) for t, v in chk]
    guardar("reto2_sensibilidad_demanda", "\n".join(md), out)
    return out



# ----------------------------------------------------------------------
# 5) Monte Carlo + reduccion de escenarios por k-means (metodo 3)
# ----------------------------------------------------------------------
MC_N, MC_SEMILLA, MC_K = 2000, 2026, (3, 6, 10, 20)


def seccion_mc():
    out, md = {}, ["# Reto 2 - metodo 3: simulacion Monte Carlo + reduccion de escenarios (k-means)\n",
                   f"El enunciado no nombra la tercera metodologia; se eligio esta. En cada bloque la intensidad eolica "
                   f"sigue una Beta ajustada por momentos (media y varianza) a los 6 niveles del paper. Se generan "
                   f"N = {MC_N} escenarios anuales (semilla {MC_SEMILLA}); por defecto las muestras de los bloques se "
                   "emparejan por rango (comonotonas), que es la dependencia implicita en los escenarios de riesgo del "
                   "paper (el w-esimo nivel de cada bloque). k-means reduce los N vectores a K escenarios con "
                   "probabilidad igual a la fraccion de muestras del grupo; cada escenario aporta 5 mercados y un grupo "
                   "de riesgo, de modo que, a diferencia del PEM, permite el CVaR del caso 4. Los planes se comparan "
                   "evaluandolos con los 30 mercados de referencia (6 niveles) y sus grupos de riesgo.\n"]
    planes_pub = {c: (parse_plan(dh.PAPER_TABLA_4[c]["ac"]), parse_plan(dh.PAPER_TABLA_4[c]["dc"])) for c in (1, 2, 3, 4)}
    # 1) distribucion
    _, _, inf6 = montecarlo_kmeans(dh.N_H, dh.CARGA, dh.VIENTO, MC_N, 6, MC_SEMILLA, True)
    md.append("## Distribucion Beta por bloque y escenarios reducidos (K = 6, comonotono)\n")
    md.append("| Bloque | Beta (a, b) | media 6 niveles | media Monte Carlo | media reducida |\n|---|---|---:|---:|---:|")
    for t in range(len(dh.N_H)):
        a, b = inf6["beta"][t]
        md.append(f"| {t + 1} | ({a:.3f}, {b:.3f}) | {inf6['media_niveles'][t]:.4f} | {inf6['media_muestral'][t]:.4f} | {inf6['media_reducida'][t]:.4f} |")
    md.append("\nProbabilidades de los 6 escenarios: " + ", ".join(f"{x:.3f}" for x in inf6["prob"]))
    out["info_K6"] = inf6
    # 2) estimador del ESW de los planes publicados (sistema del caso 3)
    S3 = dh.sistema_heydari(3)
    seg3, _ = segmentacion(dh.N_H, dh.CARGA, dh.VIENTO, True)
    md.append("\n## Exactitud del estimador del ESW (M$) de los planes publicados, sistema del caso 3\n")
    md.append("| Plan (caso) | 6 niveles (30 LP) | " + " | ".join(f"MC K={K} ({5 * K} LP)" for K in MC_K) + " |\n|---|---:|" + "---:|" * len(MC_K))
    mcc = {K: montecarlo_kmeans(dh.N_H, dh.CARGA, dh.VIENTO, MC_N, K, MC_SEMILLA, True) for K in MC_K}
    out["estimador"] = {}
    for c in (1, 2, 3, 4):
        pa, pd = planes_pub[c]
        e6 = evaluar_h(S3, pa, pd, seg3)["esw"]
        fila, celdas = [e6], []
        for K in MC_K:
            ek = evaluar_h(S3, pa, pd, mcc[K][0])["esw"]
            fila.append(ek)
            celdas.append(f"{ek:.3f} ({100 * (ek - e6) / abs(e6):+.2f} %)")
        out["estimador"][c] = fila
        md.append(f"| caso {c} | {e6:.3f} | " + " | ".join(celdas) + " |")
    # 3) plan optimo con cada K, casos 2-4
    md.append("\n## Plan optimo con cada reduccion\n")
    md.append("Arrepentimiento = objetivo (1) del plan de 6 niveles menos el del plan del metodo, ambos evaluados "
              "con los 30 mercados y los 6 grupos de riesgo de referencia (en el caso 4 incluye beta x CVaR).\n")
    md.append("| Caso | Metodo | mercados | Estado | Plan HVAC / HVDC | Inv. (M$) | Objetivo propio | CVaR propio | Objetivo con 6 niveles | Arrepentimiento (M$) | t (s) |\n|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|")
    out["opt"] = {}
    variantes = [(K, True) for K in MC_K] + [(6, False), (20, False)]
    for caso in (2, 3, 4):
        S, seg, grupos, beta = preparar(caso)
        base = resolver_h(S, seg, grupos, dh.ALFA, beta, metodo="oraculo")
        base_obj = evaluar_h(S, base.plan_ac, base.plan_dc, seg, grupos, dh.ALFA, beta)["obj"]
        md.append(f"| {caso} | 6 niveles (paper) | {len(seg)} | {base.estado} | {plan_txt(base.plan_ac)} / {plan_txt(base.plan_dc)} | "
                  f"{base.inv_total:.1f} | {base.obj:.3f} | {base.cvar if beta else 0:.3f} | {base_obj:.3f} | 0 | {base.t:.1f} |")
        out["opt"][f"{caso}-6niveles"] = dict(plan=base.plan_txt(), inv=base.inv_total, obj=base.obj, obj_ref=base_obj, t=base.t,
                                              igual_publicado=(base.plan_ac == planes_pub[caso][0] and base.plan_dc == planes_pub[caso][1]))
        for K, com in variantes:
            cm, gm, _ = montecarlo_kmeans(dh.N_H, dh.CARGA, dh.VIENTO, MC_N, K, MC_SEMILLA, com)
            r = resolver_h(S, cm, gm, dh.ALFA, beta, metodo="oraculo")
            ref_obj = evaluar_h(S, r.plan_ac, r.plan_dc, seg, grupos, dh.ALFA, beta)["obj"]
            nombre = f"MC K={K}" + ("" if com else " (bloques independientes)")
            md.append(f"| {caso} | {nombre} | {len(cm)} | {r.estado} | {plan_txt(r.plan_ac)} / {plan_txt(r.plan_dc)} | {r.inv_total:.1f} | "
                      f"{r.obj:.3f} | {r.cvar if beta else 0:.3f} | {ref_obj:.3f} | {base_obj - ref_obj:.3f} | {r.t:.1f} |")
            out["opt"][f"{caso}-K{K}{'' if com else '-indep'}"] = dict(plan=r.plan_txt(), inv=r.inv_total, obj=r.obj,
                                                                    cvar=r.cvar, obj_ref=ref_obj, regret=base_obj - ref_obj, t=r.t)
            print(f"caso {caso} {nombre}: {r.plan_txt()}  regret {base_obj - ref_obj:.3f}  t={r.t:.1f}s")
    guardar("reto2_montecarlo", "\n".join(md), out)
    return out


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if arg in ("seg", "all"):
        seccion_seg()
    if arg in ("dual", "all"):
        seccion_dual(float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0)
    if arg in ("pem", "all"):
        seccion_pem()
    if arg in ("mc", "all"):
        seccion_mc()
    if arg in ("sens", "all"):
        seccion_sens()
    print(f"\nTiempo total: {time.time() - t0:.1f} s")
