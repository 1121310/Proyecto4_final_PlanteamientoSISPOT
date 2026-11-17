# -*- coding: utf-8 -*-
"""
Datos del caso Garver de Heydari y Barforoushi (EPSR 229, 2024, 110182).

Procedencia de cada dato (nada de lo siguiente se ha ajustado para acercarse a
la Tabla 4):

  PUBLICADO EN EL PAPER
    * interes 10 %, vida util 25 anios (CRF, ec. (71)); consumo minimo 90 % de
      la oferta; costo de racionamiento = 10 x precio del 1er bloque;
      presupuesto 400 M$; coeficiente de demanda 2.2; eolicas de 500 MW en las
      barras 4 y 6 (la del 6 con una intensidad 10 % mayor); maximo 4 lineas
      HVAC y 2 HVDC por corredor; CVaR con alfa = 0.8 y beta = 5; 5 bloques de
      demanda y 6 intensidades de viento por bloque con pesos iguales (Fig. 2).
    * Tabla 3 (precios de oferta y demanda, $/MWh).
    * Caso 1: las eolicas de las barras 4 y 6 se sustituyen por unidades
      termicas de 200/150/150 MW con precios 71.5/76/80 $/MWh.
  TOMADO DE LAS FUENTES QUE EL PAPER CITA
    * Bloques de potencia de generadores y demandas "segun [1]" = Garces et al.
      2009 (datos_garver.py).
    * Datos de lineas HVAC/HVDC "segun [41]" = repositorio GP (UTP),
      http://academia.utp.edu.co/planeamiento/?p=2521 , caso "HVAC-HVDC Garver
      6-Buses, 15-lines-AC, 15-lines-DC". Costos HVAC (M$) = costos originales
      de Garver (= Tabla I de Garces / 193.08); costos HVDC (M$) transcritos de
      esa pagina (el 54.6 de la linea 2-6 coincide con la Tabla 4 de Heydari);
      la capacidad HVDC es la misma que la HVAC del corredor.
      VERIFICADO fila a fila (15 HVAC + 15 HVDC: costos, reactancias, limites,
      circuitos existentes) contra el HTML de la pagina; el snapshot esta en
      fuente_utp_garver_hvac_hvdc.csv y lo comprueba verificaciones_heydari.py.
      Unica diferencia con la pagina: ella da "Max. #L." = 5 (HVAC) y el paper
      de Heydari dice 4; se usa 4 (del paper), no es un dato de la pagina.
  INFERIDO (declararlo siempre al reportar)
    * CARGA y VIENTO: digitalizados de la Fig. 2 (el paper no los tabula; ver
      digitalizar_fig2_heydari.py).  Resolucion ~0.002 p.u.  Comprobaciones
      cruzadas con la Tabla 4 en verificaciones_heydari.py.
    * Las intensidades se ordenan de mayor a menor en cada bloque y el
      escenario omega agrupa el omega-esimo nivel de cada bloque (solo afecta
      al CVaR del caso 4).
    * "Maximo 4 lineas HVAC por corredor": se toma como 4 circuitos NUEVOS
      (variante ac_total=True: 4 contando los existentes).
    * Base de 100 MVA y barra de referencia = 1 (no afectan al resultado).
"""
import math

import datos_garver as dg
from tnep_heydari import (SistemaH, LineaAC, LineaDC, Termica, Eolica, DemandaH)

# ----------------------------------------------------------------------
# Curva de duracion y viento (Fig. 2, digitalizados)
# ----------------------------------------------------------------------
N_H = [1000, 3000, 3000, 1000, 760]                      # horas de cada bloque (suman 8760)
CARGA = [0.9143, 0.7738, 0.6350, 0.4799, 0.3345]         # p.u. de la demanda base
VIENTO = [                                               # p.u. de la capacidad instalada (barra 4)
    [0.4697, 0.3562, 0.2655, 0.1997, 0.1316, 0.0681],
    [0.4742, 0.3358, 0.2609, 0.1929, 0.1271, 0.0613],
    [0.4232, 0.2972, 0.2360, 0.1804, 0.1225, 0.0681],
    [0.5151, 0.4039, 0.3018, 0.2246, 0.1452, 0.0749],
    [0.5083, 0.4016, 0.3222, 0.2269, 0.1180, 0.0431],
]
INTERES, VIDA = 0.10, 25
CRF = INTERES * (1 + INTERES) ** VIDA / ((1 + INTERES) ** VIDA - 1)      # ec. (71) = 0.110168
ALFA, BETA = 0.8, 5.0

# ----------------------------------------------------------------------
# Lineas (repositorio GP-UTP, costos en M$ por circuito)
# ----------------------------------------------------------------------
COSTO_AC = {(1, 2): 40, (1, 3): 38, (1, 4): 60, (1, 5): 20, (1, 6): 68, (2, 3): 20, (2, 4): 40,
            (2, 5): 31, (2, 6): 30, (3, 4): 59, (3, 5): 20, (3, 6): 48, (4, 5): 63, (4, 6): 30,
            (5, 6): 61}
COSTO_DC = {(1, 2): 64.6, (1, 3): 62.6, (1, 4): 80.6, (1, 5): 44.6, (1, 6): 86.4, (2, 3): 44.6,
            (2, 4): 64.6, (2, 5): 55.6, (2, 6): 54.6, (3, 4): 80.0, (3, 5): 44.6, (3, 6): 72.6,
            (4, 5): 82.5, (4, 6): 54.6, (5, 6): 81.1}

# Tabla 3 (precios en $/MWh)
PRECIOS_GEN = {1: [66.0, 71.5, 74.8], 3: [71.5, 74.8, 77.0], 6: [77.0, 82.5, 88.0]}
PRECIOS_DEM = {1: [123.75, 90.2, 88.0], 2: [126.5, 82.5, 79.2], 3: [132.0, 88.0, 77.0],
               4: [121.0, 93.5, 82.5], 5: [123.2, 88.0, 83.6]}

# ----------------------------------------------------------------------
# Resultados publicados (Tabla 4). Plan HVAC / HVDC, M$, GWh, $/MWh
# ----------------------------------------------------------------------
PAPER_TABLA_4 = {
    1: dict(obj=442.81, ac="2-6 (4), 3-5 (3)", dc="", esw=462.64, inv=180.0, inv_ac=180.0, inv_dc=0.0,
            fosil=9627.92, eolica=0.0, consumo=9627.92, lmp=78.1, lmp_std=6.104, cvar=None),
    2: dict(obj=489.94, ac="2-6 (4), 3-5 (2), 4-6 (2), 5-6 (2)", dc="", esw=527.62, inv=342.0,
            inv_ac=342.0, inv_dc=0.0, fosil=7120.48, eolica=2213.81, consumo=9400.21, lmp=236.186,
            lmp_std=376.042, cvar=None),
    3: dict(obj=492.98, ac="2-6 (4), 3-5 (3), 3-6, 4-6 (2)", dc="2-6", esw=530.73, inv=342.6,
            inv_ac=288.0, inv_dc=54.6, fosil=7138.66, eolica=2213.81, consumo=9415.87, lmp=193.301,
            lmp_std=338.956, cvar=None),
    4: dict(obj=1652.41, ac="2-6 (4), 3-5 (3), 4-6 (3)", dc="2-6", esw=526.25, inv=324.6,
            inv_ac=270.0, inv_dc=54.6, fosil=7119.39, eolica=2210.24, consumo=9396.5, lmp=243.29,
            lmp_std=393.51, cvar=232.38),
}


def sistema_heydari(caso=2, max_ac=4, ac_total=False, max_dc=2, coef_dem=2.2, cu_factor=10.0,
                    cmax=400.0, M=1000.0, gamma=1.0e4) -> SistemaH:
    """Garver de Heydari. caso 1: solo termicas (eolicas sustituidas por termicas);
    caso 2: con eolica, solo HVAC; casos 3 y 4: con eolica y candidatas HVDC."""
    ac = []
    for c, (o, r, x, f, _) in enumerate(dg._TABLA_I):
        nex = 1 if (o, r) in dg._EXISTENTES else 0
        ncand = (max_ac - nex) if ac_total else max_ac
        for t in range(nex + ncand):
            ac.append(LineaAC(o, r, x, float(f), float(COSTO_AC[(o, r)]), existente=(t < nex),
                              corredor=c, copia=t))
    dc = []
    if caso >= 3:
        for c, (o, r, x, f, _) in enumerate(dg._TABLA_I):
            for t in range(max_dc):
                dc.append(LineaDC(o, r, float(f), float(COSTO_DC[(o, r)]), corredor=c, copia=t))
    term = [Termica(b, list(zip(dg._GENS[b][0], PRECIOS_GEN[b]))) for b in (1, 3, 6)]
    eol = []
    if caso == 1:
        for b in (4, 6):
            term.append(Termica(b, [(200.0, 71.5), (150.0, 76.0), (150.0, 80.0)]))
    else:
        eol = [Eolica(4, 500.0, 1.0), Eolica(6, 500.0, 1.1)]
    dems = [DemandaH(b, [(coef_dem * mw, p) for mw, p in zip(dg._DEMS[b][0], PRECIOS_DEM[b])],
                     cu_factor=cu_factor) for b in sorted(dg._DEMS)]
    return SistemaH("Garver-Heydari", [1, 2, 3, 4, 5, 6], slack=1, ac=ac, dc=dc, term=term, eol=eol,
                    dems=dems, cmax=cmax, crf=CRF, M=M, Gamma=gamma, GammaDC=gamma)
