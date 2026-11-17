# -*- coding: utf-8 -*-
"""
Paso 1: verificar la transcripcion de datos recalculando, a partir de las
tablas transcritas, los costos de inversion de TODOS los planes publicados.
Si una cifra de una tabla fuera mal transcrita, alguna de estas igualdades
fallaria.
"""
import datos_garver as G
import datos_rts24 as R
from sistema import parse_plan

TOL = 0.006   # los papers redondean a 2 decimales (M EUR)


def chequea(sistema, nombre, plan_txt, inv_publicada):
    plan = parse_plan(plan_txt)
    calc = sistema.costo_plan(plan)
    ok = abs(calc - inv_publicada) <= TOL
    print(f"  {'OK ' if ok else 'ERR'} {nombre:<28s} plan = {plan_txt:<42s} "
          f"calculado = {calc:8.3f}  publicado = {inv_publicada:6.2f}")
    assert ok, f"costo de plan no coincide: {nombre}"


def main():
    print("== Garver: planes de la Tabla IV y V ==")
    g = G.sistema_garver()
    for k, v in G.PAPER_TABLA_IV.items():
        chequea(g, k, v["plan"], v["inv"])
    for k, v in G.PAPER_TABLA_V.items():
        chequea(g, f"Contingencia {k[0]}-{k[1]}", v["plan"], v["inv"])

    print("== RTS-24: planes de las Tablas VIII, X y Cuadro 6.9 ==")
    r = R.sistema_rts24(bloques="impresos")
    for k, v in R.PAPER_TABLA_VIII.items():
        chequea(r, k, v["plan"], v["inv"])
    for k, v in R.PAPER_CUADRO_6_9.items():
        chequea(r, f"Contingencia {k[0]}-{k[1]}", v["plan"], v["inv"])
    chequea(r, "Multicontingencia (Tabla X)", R.PAPER_TABLA_X["plan"], R.PAPER_TABLA_X["inv"])

    print("== Consistencia de las demas tablas ==")
    # Garver: cifras de la Fig. 1
    gen_bus = {}
    for ge in g.gens:
        gen_bus[ge.bus] = sum(t for t, _ in ge.bloques)
    dem_bus = {d.bus: sum(t for t, _ in d.bloques) for d in g.dems}
    assert gen_bus == {1: 150, 3: 360, 6: 600}, gen_bus
    assert dem_bus == {1: 80, 2: 240, 3: 40, 4: 160, 5: 240}, dem_bus
    print("  OK  Garver: generacion 150/360/600 MW y demandas 80/240/40/160/240 MW (Fig. 1)")
    assert len(g.lineas) == 15 * 3 and len(g.candidatas) == 39
    assert sum(1 for L in g.lineas if L.existente) == 6
    print(f"  OK  Garver: 15 corredores, 6 circuitos existentes, {len(g.candidatas)} candidatos")

    # RTS: Cuadro A.5 (demandas maximas) y consistencia con 3x RTS original
    dem_rts = {d.bus: sum(t for t, _ in d.bloques) for d in r.dems}
    assert dem_rts == R.CUADRO_A5, (dem_rts, R.CUADRO_A5)
    print(f"  OK  RTS-24: demandas por barra = Cuadro A.5; total = {sum(dem_rts.values())} MW (3 x 2850)")
    gen_rts = {ge.bus: sum(t for t, _ in ge.bloques) for ge in r.gens}
    esperado = {1: 576, 2: 576, 7: 900, 13: 1773, 15: 645, 16: 465, 18: 1200,
                21: 1200, 22: 900, 23: 1980}     # 3 x (192,192,300,591,215,155,400,400,300,660)
    assert gen_rts == esperado, gen_rts
    print(f"  OK  RTS-24: generacion total = {sum(gen_rts.values())} MW (3 x 3405), por barra = 3 x RTS original")
    # bloques exactos (80/10/10 y 40/30/30 %) vs impresos: mismos totales, |dif| <= 0.6 MW
    rx = R.sistema_rts24(bloques="exactos")
    dmax = 0.0
    for a, b in list(zip(r.dems, rx.dems)) + list(zip(r.gens, rx.gens)):
        assert a.bus == b.bus
        assert abs(sum(t for t, _ in a.bloques) - sum(t for t, _ in b.bloques)) < 1e-9
        assert [pr for _, pr in a.bloques] == [pr for _, pr in b.bloques]
        dmax = max(dmax, max(abs(t1 - t2) for (t1, _), (t2, _) in zip(a.bloques, b.bloques)))
    assert dmax <= 0.6 + 1e-9, dmax
    print(f"  OK  RTS-24: bloques 'exactos' vs 'impresos': mismos totales y precios; desviacion maxima = {dmax:.2f} MW (redondeo)")
    n_cor = len({(L.o, L.r) for L in r.lineas})
    n_exist_cor = len({(L.o, L.r) for L in r.lineas if L.existente})
    n_exist_circ = sum(1 for L in r.lineas if L.existente)
    assert (n_cor, n_exist_cor, n_exist_circ) == (41, 34, 38), (n_cor, n_exist_cor, n_exist_circ)
    print(f"  OK  RTS-24: {n_cor} corredores (34 existentes con {n_exist_circ} circuitos + 7 nuevos), "
          f"{len(r.candidatas)} circuitos candidatos")

    # Pesos de escenarios
    for nombre, esc in [("Garver base", G.escenarios_base()),
                        ("Garver cont.", G.escenarios_contingencia((1, 5))),
                        ("RTS sin cont.", R.escenarios_sin_contingencia()),
                        ("RTS cont. simple", R.escenarios_contingencia_simple((11, 13)))]:
        s = sum(e.peso for e in esc)
        assert abs(s - 1.0) < 1e-12, (nombre, s)
        print(f"  OK  {nombre}: {len(esc)} escenarios, suma de pesos = {s:.4f}")
    esc = R.escenarios_multicontingencia()
    s = sum(e.peso for e in esc)
    print(f"  !!  RTS multicontingencia: {len(esc)} escenarios, suma de pesos = {s:.4f} "
          f"(0.25 + 21 x 0.0357, tal como se imprime en la Tabla IX; NO se renormaliza)")
    print("\nTodas las verificaciones de datos pasaron.")


if __name__ == "__main__":
    main()
