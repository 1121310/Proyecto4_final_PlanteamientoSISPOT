# -*- coding: utf-8 -*-
"""
IEEE RTS de 24 barras (caso de estudio del paper, Seccion IV).

Fuentes (todas incluidas en la carpeta del proyecto):
  * Lineas (X, capacidad, costo, n. de circuitos existentes): tesis de Garces,
    Apendice A.3, Cuadro A.3 (34 corredores existentes + 7 nuevos).
  * Generadores y demandas (bloques y precios): paper Tabla VI = tesis Cuadro 6.6.
  * Escenarios: paper Tablas VII y IX = tesis Cuadros 6.7 y 6.10.
"""
from sistema import (Sistema, Generador, Demanda, Escenario, construir_lineas)

# (o, r, capacidad MW, X p.u., costo 10^3 EUR, n. circuitos existentes)
_CUADRO_A3 = [
    (1, 2, 175, 0.0139, 337.89, 1),
    (1, 3, 175, 0.2112, 6194.65, 1),
    (1, 5, 175, 0.0845, 2477.86, 1),
    (2, 4, 175, 0.1267, 3716.79, 1),
    (2, 6, 175, 0.1920, 5631.50, 1),
    (3, 9, 175, 0.1190, 3491.53, 1),
    (3, 24, 400, 0.0839, 9654.00, 1),
    (4, 9, 175, 0.1037, 3041.01, 1),
    (5, 10, 175, 0.0883, 2590.49, 1),
    (6, 10, 175, 0.0605, 1802.08, 1),
    (7, 8, 175, 0.0614, 1802.08, 1),
    (8, 9, 175, 0.1651, 4843.09, 1),
    (8, 10, 175, 0.1651, 4843.09, 1),
    (9, 11, 400, 0.0839, 9654.00, 1),
    (9, 12, 400, 0.0839, 9654.00, 1),
    (10, 11, 400, 0.0839, 9654.00, 1),
    (10, 12, 400, 0.0839, 9654.00, 1),
    (11, 13, 500, 0.0476, 6371.64, 1),
    (11, 14, 500, 0.0418, 5599.32, 1),
    (12, 13, 500, 0.0476, 6371.64, 1),
    (12, 23, 500, 0.0966, 12936.36, 1),
    (13, 23, 500, 0.0865, 11584.80, 1),
    (14, 16, 500, 0.0389, 5213.16, 1),
    (15, 16, 500, 0.0173, 2316.96, 1),
    (15, 21, 500, 0.0490, 6564.72, 2),
    (15, 24, 500, 0.0519, 6950.88, 1),
    (16, 17, 500, 0.0259, 3475.44, 1),
    (16, 19, 500, 0.0231, 3089.28, 1),
    (17, 18, 500, 0.0144, 1930.80, 1),
    (17, 22, 500, 0.1053, 14094.84, 1),
    (18, 21, 500, 0.0259, 3475.44, 2),
    (19, 20, 500, 0.0396, 5309.70, 2),
    (20, 23, 500, 0.0216, 2896.20, 2),
    (21, 22, 500, 0.0678, 9074.76, 1),
    # corredores nuevos (0 circuitos existentes)
    (1, 8, 175, 0.1344, 2450.00, 0),
    (2, 8, 175, 0.1267, 2310.00, 0),
    (6, 7, 175, 0.1920, 3500.00, 0),
    (13, 14, 500, 0.0447, 3720.00, 0),
    (14, 23, 500, 0.0620, 5160.00, 0),
    (16, 23, 500, 0.0822, 6840.00, 0),
    (19, 23, 500, 0.0606, 5040.00, 0),
]

# Tabla VI: bus -> (tamanos MW, precios EUR/MWh)
_GENS = {
    1:  ([230, 173, 173],     [95, 97, 98]),
    2:  ([230, 173, 173],     [96, 99, 100]),
    7:  ([360, 270, 270],     [96, 98, 100]),
    13: ([709, 532, 532],     [80, 82, 85]),
    15: ([258, 193.5, 193.5], [82, 84, 86]),
    16: ([186, 139.5, 139.5], [77, 79, 80]),
    18: ([480, 360, 360],     [73, 76, 78]),
    21: ([480, 360, 360],     [74, 76, 78]),
    22: ([360, 270, 270],     [79, 81, 83]),
    23: ([792, 594, 594],     [78, 80, 82]),
}
_DEMS = {
    1:  ([259, 32.5, 32.5], [99, 80, 78]),
    2:  ([233, 29, 29],     [98, 81, 77]),
    3:  ([432, 54, 54],     [100, 85, 83]),
    4:  ([178, 22, 22],     [99, 86, 83]),
    5:  ([171, 21, 21],     [100, 85, 82]),
    6:  ([326, 41, 41],     [99, 85, 81]),
    7:  ([300, 37.5, 37.5], [100, 86, 83]),
    8:  ([411, 51, 51],     [93, 85, 82]),
    9:  ([420, 52.5, 52.5], [99, 82, 81]),
    10: ([468, 58.5, 58.5], [100, 80, 78]),
    13: ([636, 79.5, 79.5], [92, 79, 77]),
    14: ([466, 58, 58],     [90, 78, 76]),
    15: ([761, 95, 95],     [87, 77, 76]),
    16: ([240, 30, 30],     [84, 75, 73]),
    18: ([799, 100, 100],   [91, 76, 72]),
    19: ([435, 54, 54],     [94, 75, 74]),
    20: ([307, 38.5, 38.5], [95, 80, 77]),
}

# Demandas maximas por barra (tesis, Cuadro A.5) - solo para verificar la transcripcion
CUADRO_A5 = {1: 324, 2: 291, 3: 540, 4: 222, 5: 213, 6: 408, 7: 375, 8: 513,
             9: 525, 10: 585, 13: 795, 14: 582, 15: 951, 16: 300, 18: 999,
             19: 543, 20: 384}

TODAS = list(range(1, 25))
SO = [1, 3, 4, 5, 9]          # suroeste
SE = [2, 6, 7, 8, 10]         # sureste
NO = [15, 16, 18]             # noroeste
NE = [13, 14, 19, 20]         # noreste
SUR = list(range(1, 11))      # barras 1-10
NORTE = list(range(11, 25))   # barras 11-24


def _f(*grupos):
    d = {}
    for valor, barras in grupos:
        for b in barras:
            d[b] = valor
    return d


def escenarios_sin_contingencia():
    """Tabla VII."""
    return [
        Escenario("1", 0.50, _f((1.0, TODAS))),
        Escenario("2", 0.25, _f((1.1, SUR), (1.0, NORTE))),
        Escenario("3", 0.25, _f((1.0, SUR), (1.2, NORTE))),
    ]


def escenarios_contingencia_simple(corredor):
    """Cuadro 6.9 de la tesis: 3 escenarios originales (pesos 0.40/0.25/0.25) +
    un escenario de peso 0.10, factores 1 y la salida de una linea."""
    sc = escenarios_sin_contingencia()
    sc[0] = Escenario("1", 0.40, dict(sc[0].factores))
    sc.append(Escenario(f"4 (cont. {corredor[0]}-{corredor[1]})", 0.10,
                        _f((1.0, TODAS)), linea_fuera=corredor))
    return sc


def escenarios_multicontingencia():
    """Tabla IX (22 escenarios). Pesos tal como se imprimen (suman 0.9997)."""
    w = 0.0357
    s = []
    s.append(Escenario("1", 0.25, _f((1.0, TODAS))))
    s.append(Escenario("2", w, _f((1.1, SO), (1.0, SE + NORTE))))
    s.append(Escenario("3", w, _f((1.0, SO + NORTE), (1.1, SE))))
    s.append(Escenario("4", w, _f((1.0, SUR + NO), (1.2, NE))))
    s.append(Escenario("5", w, _f((1.0, SUR + NE), (1.2, NO))))
    s.append(Escenario("6", w, _f((1.1, SUR), (1.0, NORTE))))
    s.append(Escenario("7", w, _f((1.1, SO), (1.0, SE + NO), (1.2, NE))))
    s.append(Escenario("8", w, _f((1.1, SO), (1.0, SE + NE), (1.2, NO))))
    s.append(Escenario("9", w, _f((1.0, SO + NO), (1.1, SE), (1.2, NE))))
    s.append(Escenario("10", w, _f((1.0, SO + NE), (1.1, SE), (1.2, NO))))
    s.append(Escenario("11", w, _f((1.0, SUR), (1.2, NORTE))))
    s.append(Escenario("12", w, _f((1.1, SUR), (1.2, NE), (1.0, NO))))
    s.append(Escenario("13", w, _f((1.1, SUR), (1.0, NE), (1.2, NO))))
    s.append(Escenario("14", w, _f((1.1, SO), (1.0, SE), (1.2, NORTE))))
    s.append(Escenario("15", w, _f((1.0, SO), (1.1, SE), (1.2, NORTE))))
    s.append(Escenario("16", w, _f((1.1, SUR), (1.2, NORTE))))
    for n, cor in zip(range(17, 23),
                      [(11, 13), (11, 14), (12, 13), (12, 23), (14, 16), (15, 24)]):
        s.append(Escenario(str(n), w, _f((1.0, TODAS)), linea_fuera=cor))
    return s


# Fracciones de la capacidad/demanda total que se reparten en los 3 bloques.
# Las tablas impresas (Tabla VI) son una version REDONDEADA de estas fracciones
# (p. ej. barra 1: 0.8*324 = 259.2 -> "259", 0.1*324 = 32.4 -> "32.5"; 0.4*576 = 230.4 -> "230"),
# con ajustes para conservar la suma. En Garver los bloques son exactos (72/4/4 = 90/5/5 %).
FRAC_DEM = (0.8, 0.1, 0.1)
FRAC_GEN = (0.4, 0.3, 0.3)


def sistema_rts24(cmax=23.0, bloques="exactos") -> Sistema:
    """bloques='impresos': tamanos de bloque tal como aparecen en la Tabla VI.
       bloques='exactos' : mismos totales por barra, repartidos en 80/10/10 % (demanda)
                           y 40/30/30 % (generacion). Es una INFERENCIA, no un dato del
                           paper; se justifica porque reproduce las 4 cifras de bienestar
                           medio de la Tabla VIII con error <= 0.003 M EUR (los datos
                           impresos las desplazan +0.19 M EUR)."""
    lineas = construir_lineas([(o, r, x, f, c, n) for (o, r, f, x, c, n) in _CUADRO_A3],
                              max_lineas=3)
    gens, dems = [], []
    for b, (t, p) in _GENS.items():
        tam = [f * sum(t) for f in FRAC_GEN] if bloques == "exactos" else t
        gens.append(Generador(b, list(zip(tam, p))))
    for b, (t, p) in _DEMS.items():
        tam = [f * sum(t) for f in FRAC_DEM] if bloques == "exactos" else t
        dems.append(Demanda(b, list(zip(tam, p))))
    return Sistema("IEEE-RTS-24" + ("" if bloques == "exactos" else " (bloques impresos)"),
                   list(range(1, 25)), slack=1, lineas=lineas, gens=gens, dems=dems, cmax=cmax)


# ---------------- Resultados publicados ----------------
PAPER_TABLA_VIII = {
    "Bilevel":           dict(plan="6-10, 10-11, 11-13, 14-23",     inv=22.99, sw=None,   swm=742.59,   cpu=89.51),
    "Min-costos Esc. 1": dict(plan="6-7, 7-8, 13-14",               inv=9.02,  sw=587.83, swm=-2430.44, cpu=2.68),
    "Min-costos Esc. 2": dict(plan="6-7 (2), 7-8, 2-8, 11-13",      inv=17.48, sw=606.54, swm=687.17,   cpu=40.29),
    "Min-costos Esc. 3": dict(plan="6-10, 7-8 (2), 13-14",          inv=9.13,  sw=621.62, swm=-1399.75, cpu=2.68),
}

PAPER_CUADRO_6_9 = {   # contingencia simple, tesis (CPU en minutos)
    (11, 13): dict(plan="6-10, 10-11, 11-13, 14-23",      inv=22.99, swm=740.39, cpu_min=4.59),
    (11, 14): dict(plan="6-10, 10-11, 11-13, 14-23",      inv=22.99, swm=740.24, cpu_min=4.20),
    (12, 13): dict(plan="6-10, 10-11, 11-13, 14-23",      inv=22.99, swm=740.00, cpu_min=4.96),
    (12, 23): dict(plan="5-10, 6-10, 7-8 (2), 10-11, 14-23", inv=22.81, swm=713.02, cpu_min=16.19),
    (14, 16): dict(plan="6-10, 7-8, 10-12, 14-16, 13-14", inv=22.19, swm=725.42, cpu_min=6.21),
    (15, 24): dict(plan="3-9, 6-10, 7-8 (2), 10-12, 13-14", inv=22.27, swm=711.89, cpu_min=7.38),
}

PAPER_TABLA_X = dict(plan="3-9, 6-10, 7-8 (2), 10-11, 11-13, 13-14", inv=28.64, swm=717.97, cpu_h=5.43)
