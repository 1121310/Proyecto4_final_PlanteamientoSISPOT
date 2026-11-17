# -*- coding: utf-8 -*-
"""
Paso 2: el nivel inferior como LP primal y como LP dual independientes.
Con el plan x FIJO, el optimo primal (5)-(19) debe igualar el optimo del dual
(21)-(34) (dualidad fuerte). Se prueban varios planes y escenarios, incluidos
casos con racionamiento (precios ~ 100 x oferta), donde la cota Gamma_max=1000
sobre phi podria ser activa.
"""
import datos_garver as G
import datos_rts24 as R
from sistema import parse_plan
import tnep

TOL_REL = 1e-7


def prueba(S, nombre_sis, planes, escenarios):
    print(f"\n=== {nombre_sis} ===")
    print(f"{'plan':<36s}{'esc':>4s}{'racion[MW]':>12s}{'SW primal':>16s}{'SW dual':>16s}"
          f"{'dif.rel':>10s}{'max|phi|':>10s}{'dual c/Gamma':>14s}")
    n_ok = 0
    for txt in planes:
        plan = parse_plan(txt) if txt else {}
        sw_prom, sws, racs = tnep.evaluar_plan(S, plan, escenarios)
        for w, sc in enumerate(escenarios):
            dobj, maxphi = tnep.dual_lp(S, plan, sc, gamma=None)
            dg, _ = tnep.dual_lp(S, plan, sc, gamma=S.Gamma)
            rel = abs(sws[w] - dobj) / max(1.0, abs(sws[w]))
            assert rel < TOL_REL, f"dualidad fuerte falla: {txt} esc {sc.nombre}: {sws[w]} vs {dobj}"
            n_ok += 1
            gtxt = "factible=" + ("= SW" if (dg is not None and abs(dg - dobj) < 1e-4 * max(1, abs(dobj)))
                                  else ("distinto" if dg is not None else "INFACTIBLE"))
            print(f"{(txt or '(sin expansion)'):<36s}{sc.nombre:>4s}{racs[w]:12.2f}{sws[w]:16.4f}"
                  f"{dobj:16.4f}{rel:10.1e}{maxphi:10.1f}{gtxt:>14s}")
    print(f"--> {n_ok} pruebas primal=dual OK")


def main():
    g = G.sistema_garver()
    esc = G.escenarios_base()
    prueba(g, "Garver", ["", "3-5, 2-3, 4-6 (3)", "3-5, 2-3, 2-6 (2)", "3-5, 2-6, 4-6 (2)",
                         "3-5, 4-6 (3)", "1-6, 2-6, 3-6"], esc)
    r = R.sistema_rts24()
    esc_r = R.escenarios_sin_contingencia()
    prueba(r, "RTS-24", ["", "6-10, 10-11, 11-13, 14-23", "6-7, 7-8, 13-14", "6-7 (2), 7-8, 2-8, 11-13"],
           esc_r)
    # un escenario con contingencia
    prueba(r, "RTS-24 contingencia 14-16", ["6-10, 10-11, 11-13, 14-23"],
           R.escenarios_contingencia_simple((14, 16))[3:])
    print("\nLas pruebas de dualidad primal = dual pasaron.")


if __name__ == "__main__":
    main()
