# -*- coding: utf-8 -*-
"""
Sistema de Garver de 6 barras - datos transcritos del paper (Tablas I, II, III)
y de la tesis de Garces (Cuadros 6.1, 6.2, 6.3, identicos).
"""
from sistema import (Sistema, Generador, Demanda, Escenario, construir_lineas)

# (o, r, X [p.u.], capacidad [MW], inversion [10^3 EUR], n. circuitos existentes)
# Existentes (6): 1-2, 1-4, 1-5, 2-3, 2-4, 3-5  (Fig. 1 del paper, "six existing lines")
_EXISTENTES = {(1, 2), (1, 4), (1, 5), (2, 3), (2, 4), (3, 5)}
_TABLA_I = [
    (1, 2, 0.40, 100, 7723.20),
    (1, 3, 0.38, 100, 7337.04),
    (1, 4, 0.60, 80, 11584.80),
    (1, 5, 0.20, 100, 3861.60),
    (1, 6, 0.68, 70, 13129.44),
    (2, 3, 0.20, 100, 3861.60),
    (2, 4, 0.40, 100, 7723.20),
    (2, 5, 0.31, 100, 5985.48),
    (2, 6, 0.30, 100, 5792.40),
    (3, 4, 0.59, 82, 11391.72),
    (3, 5, 0.20, 100, 3861.60),
    (3, 6, 0.48, 100, 9267.84),
    (4, 5, 0.63, 75, 12164.04),
    (4, 6, 0.30, 100, 5792.40),
    (5, 6, 0.61, 78, 11777.88),
]

# Tabla II: bloques (tamano MW) y precios (EUR/MWh)
_GENS = {
    1: ([60, 45, 45], [60, 65, 68]),
    3: ([144, 108, 108], [65, 68, 70]),
    6: ([240, 180, 180], [70, 75, 80]),
}
_DEMS = {
    1: ([72, 4, 4], [112.5, 82, 80]),
    2: ([216, 12, 12], [115, 75, 72]),
    3: ([36, 2, 2], [120, 80, 70]),
    4: ([144, 8, 8], [110, 85, 75]),
    5: ([216, 12, 12], [112, 80, 76]),
}


def _factores(*grupos):
    f = {}
    for valor, barras in grupos:
        for b in barras:
            f[b] = valor
    return f


def escenarios_base():
    """Tabla III."""
    return [
        Escenario("1", 0.50, _factores((1.0, [1, 2, 3, 4, 5]))),
        Escenario("2", 0.25, _factores((1.0, [1, 3, 4]), (1.2, [2, 5]))),
        Escenario("3", 0.25, _factores((1.0, [2, 3, 5]), (1.4, [1, 4]))),
    ]


def escenarios_contingencia(corredor):
    """Seccion III-B: escenario adicional de peso 0.10, factores = 1 y la
    salida de una linea; el peso del escenario 1 baja a 0.40."""
    sc = escenarios_base()
    sc[0] = Escenario("1", 0.40, dict(sc[0].factores))
    sc.append(Escenario(f"4 (cont. {corredor[0]}-{corredor[1]})", 0.10,
                        _factores((1.0, [1, 2, 3, 4, 5])), linea_fuera=corredor))
    return sc


def sistema_garver(cmax=30.0) -> Sistema:
    tabla = [(o, r, x, f, c, 1 if (o, r) in _EXISTENTES else 0)
             for (o, r, x, f, c) in _TABLA_I]
    lineas = construir_lineas(tabla, max_lineas=3)
    gens = [Generador(b, list(zip(t, p))) for b, (t, p) in _GENS.items()]
    dems = [Demanda(b, list(zip(t, p))) for b, (t, p) in _DEMS.items()]
    return Sistema("Garver-6", [1, 2, 3, 4, 5, 6], slack=1, lineas=lineas,
                   gens=gens, dems=dems, cmax=cmax)


# Resultados publicados (paper Tabla IV, V; tesis Cuadros 6.4, 6.5)
PAPER_TABLA_IV = {
    "Bilevel":            dict(plan="3-5, 2-3, 4-6 (3)", inv=25.10, sw=None,   swm=292.71,   cpu=6.70),
    "Min-costos Esc. 1":  dict(plan="3-5, 2-3, 2-6 (2)", inv=19.31, sw=275.60, swm=-2613.02, cpu=0.68),
    "Min-costos Esc. 2":  dict(plan="3-5, 2-3, 4-6 (3)", inv=25.10, sw=291.50, swm=292.71,   cpu=0.70),
    "Min-costos Esc. 3":  dict(plan="3-5, 2-6, 4-6 (2)", inv=21.24, sw=286.62, swm=-423.57,  cpu=0.52),
    # Solo en la tesis (Cuadro 6.4):
    "Basica y centralizada (tesis)": dict(plan="3-5, 4-6 (3)", inv=21.24, sw=None, swm=-507.84, cpu=None),
}

PAPER_TABLA_V = {   # contingencia simple (paper Tabla V)
    (1, 2): dict(plan="2-3, 3-5, 4-6 (3)",         inv=25.10, swm=292.71, cpu=4.19),
    (1, 4): dict(plan="2-3, 3-5, 4-6 (3)",         inv=25.10, swm=292.71, cpu=1.65),
    (1, 5): dict(plan="1-5, 2-3, 3-5, 4-6 (3)",    inv=28.96, swm=292.72, cpu=5.09),
    (2, 3): dict(plan="2-3, 3-5, 4-6 (3)",         inv=25.10, swm=292.53, cpu=6.05),
    (2, 4): dict(plan="2-3, 3-5, 4-6 (3)",         inv=25.10, swm=292.71, cpu=8.80),
    (3, 5): dict(plan="2-3, 3-5 (2), 4-6 (3)",     inv=28.96, swm=292.70, cpu=23.30),
}
