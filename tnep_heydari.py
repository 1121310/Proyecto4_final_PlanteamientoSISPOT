# -*- coding: utf-8 -*-
"""
TNEP bilevel con lineas HVAC/HVDC, generacion eolica y riesgo (CVaR):

    R. Heydari, T. Barforoushi, "HVDC/HVAC transmission network expansion
    planning in electricity markets in the presence of wind resources",
    Electric Power Systems Research 229 (2024) 110182.
    (Ecuaciones (1)-(71) referenciadas en los comentarios.)

Es la extension del modelo de Garces et al. (2009) de tnep.py:
  * unidades eolicas con costo 0 y disponibilidad por caso, ec. (19);
  * lineas HVDC con flujo controlable |f| <= x f_max (sin angulo), ec. (14)-(16),
    cuyo dual es (35)-(36) con la linealizacion (64)-(66);
  * CVaR en el nivel superior, ec. (1)-(2), (8);
  * costo de racionamiento = cu_factor x precio del primer bloque (10 en Heydari).

Estructura de la incertidumbre (independiente del metodo de modelado):
  Caso  : un problema de mercado (bloque de la curva de duracion, nivel de
          carga, intensidad eolica) con `horas` = horas/ano ponderadas por su
          probabilidad  ->  E[bienestar anual] = sum_c horas_c * SW_c.
  Grupo : escenario de riesgo (una trayectoria anual completa) con
          probabilidad delta y los casos (con sus horas de bloque) que la
          componen; solo se usa si beta > 0 (CVaR).

Unidades: potencias MW, precios $/MWh, inversion M$ (costo TOTAL por circuito),
bienestar y objetivo en M$/ano.
"""
import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import gurobipy as gp
from gurobipy import GRB

from tnep import _nuevo_modelo, _estado

HORAS_ANO = 8760.0


# ----------------------------------------------------------------------
# datos
# ----------------------------------------------------------------------
@dataclass
class LineaAC:
    o: int
    r: int
    x: float
    fmax: float
    costo: float          # M$ por circuito
    existente: bool
    corredor: int
    copia: int


@dataclass
class LineaDC:
    o: int
    r: int
    fmax: float
    costo: float          # M$ por circuito
    corredor: int
    copia: int


@dataclass
class Termica:
    bus: int
    bloques: List[Tuple[float, float]]      # (MW, $/MWh)


@dataclass
class Eolica:
    bus: int
    cap: float                              # MW instalados
    factor: float = 1.0                     # intensidad relativa al bus de referencia


@dataclass
class DemandaH:
    bus: int
    bloques: List[Tuple[float, float]]      # (MW a carga 1 p.u., $/MWh)
    cu_factor: float = 10.0                 # Heydari: 10 x precio del 1er bloque

    @property
    def cU(self) -> float:
        return self.cu_factor * self.bloques[0][1]


@dataclass
class SistemaH:
    nombre: str
    barras: List[int]
    slack: int
    ac: List[LineaAC]
    dc: List[LineaDC]
    term: List[Termica]
    eol: List[Eolica]
    dems: List[DemandaH]
    cmax: float                    # presupuesto sobre el costo TOTAL (M$)
    crf: float
    Sbase: float = 100.0
    dmin_frac: float = 0.90
    M: float = 1000.0
    Gamma: float = 1.0e4           # cota de phi (Heydari no la indica)
    GammaDC: float = 1.0e4

    @property
    def cand_ac(self):
        return [k for k, L in enumerate(self.ac) if not L.existente]

    @property
    def cand_dc(self):
        return list(range(len(self.dc)))


@dataclass
class Caso:
    nombre: str
    bloque: int
    horas: float        # horas/ano ponderadas por la probabilidad del caso
    carga: float        # nivel de carga (p.u. de la demanda base)
    viento: float       # intensidad eolica (p.u. de la capacidad instalada)


@dataclass
class Grupo:
    nombre: str
    delta: float
    miembros: List[Tuple[int, float]]       # (indice de caso, horas completas del bloque)


def costo_plan(S: SistemaH, plan_ac: dict, plan_dc: dict) -> Tuple[float, float, float]:
    """(total, HVAC, HVDC) en M$ de un plan {(o,r): n_circuitos_nuevos}."""
    ca = cd = 0.0
    for (o, r), n in (plan_ac or {}).items():
        c = [L.costo for L in S.ac if not L.existente and (L.o, L.r) == (o, r)]
        if not c:
            raise ValueError(f"corredor HVAC {(o, r)} sin candidatas")
        ca += n * c[0]
    for (o, r), n in (plan_dc or {}).items():
        c = [L.costo for L in S.dc if (L.o, L.r) == (o, r)]
        if not c:
            raise ValueError(f"corredor HVDC {(o, r)} sin candidatas")
        cd += n * c[0]
    return ca + cd, ca, cd


def plan_txt(plan: dict) -> str:
    return ", ".join(f"{o}-{r}" + (f" ({n})" if n > 1 else "") for (o, r), n in sorted(plan.items())) or "-"


def parse_plan(texto: str) -> Dict[Tuple[int, int], int]:
    """'2-6 (4), 3-5 (3)' o '2–6 (4), 3–5 (3)' -> {(2,6):4,(3,5):3}; '' / '-' -> {}"""
    plan = {}
    for tok in texto.replace("–", "-").split(","):
        tok = tok.strip()
        if not tok or tok == "-":
            continue
        n = 1
        if "(" in tok:
            tok, rest = tok.split("(")
            n = int(rest.strip(") "))
        a, b = tok.strip().split("-")
        plan[(int(a), int(b))] = n
    return plan


def _x_fijo(lineas, plan, es_ac):
    x, cont = {}, {}
    for k, L in enumerate(lineas):
        if es_ac and L.existente:
            x[k] = 1.0
            continue
        c = cont.get((L.o, L.r), 0)
        x[k] = 1.0 if c < plan.get((L.o, L.r), 0) else 0.0
        cont[(L.o, L.r)] = c + 1
    for (o, r), n in plan.items():
        cand = sum(1 for L in lineas if (L.o, L.r) == (o, r) and (not es_ac or not L.existente))
        if n > cand:
            raise ValueError(f"plan pide {n} circuitos en {o}-{r}; maximo {cand}")
    return x


def _a_plan(lineas, xv, es_ac):
    plan = {}
    for k, L in enumerate(lineas):
        if es_ac and L.existente:
            continue
        if xv[k] > 0.5:
            plan[(L.o, L.r)] = plan.get((L.o, L.r), 0) + 1
    return plan


# ----------------------------------------------------------------------
# Nivel inferior: mercado de un caso, ec. (9)-(28) con (59)-(63)
# ----------------------------------------------------------------------
def _idx_bus(S: SistemaH):
    gb = {s: [] for s in S.barras}
    eb = {s: [] for s in S.barras}
    db = {s: [] for s in S.barras}
    for i, g in enumerate(S.term):
        gb[g.bus].append(i)
    for y, e in enumerate(S.eol):
        eb[e.bus].append(y)
    for j, d in enumerate(S.dems):
        db[d.bus].append(j)
    return gb, eb, db


def _dmax_blq(S, cs: Caso, j, h):
    return cs.carga * S.dems[j].bloques[h][0]


def _dmax_tot(S, cs: Caso, j):
    return sum(_dmax_blq(S, cs, j, h) for h in range(len(S.dems[j].bloques)))


def _pw(S, cs: Caso, y):
    e = S.eol[y]
    return e.cap * e.factor * cs.viento


def add_primal_h(m, S: SistemaH, cs: Caso, xac: dict, xdc: dict, tag: str, M=None):
    """x = float (plan fijo: red exacta) o gp.Var (MILP, linealizacion big-M)."""
    M = S.M if M is None else M
    gb, eb, db = _idx_bus(S)
    g = {(i, b): m.addVar(0.0, G.bloques[b][0], name=f"g_{tag}_{i}_{b}")
         for i, G in enumerate(S.term) for b in range(len(G.bloques))}                    # (17)
    gw = {y: m.addVar(0.0, _pw(S, cs, y), name=f"gw_{tag}_{y}") for y in range(len(S.eol))}   # (19)
    d = {(j, h): m.addVar(0.0, _dmax_blq(S, cs, j, h), name=f"d_{tag}_{j}_{h}")
         for j, D in enumerate(S.dems) for h in range(len(D.bloques))}                    # (18)
    r = {j: m.addVar(0.0, _dmax_tot(S, cs, j), name=f"r_{tag}_{j}") for j in range(len(S.dems))}   # (20)
    th = {s: m.addVar(-math.pi, math.pi, name=f"th_{tag}_{s}") for s in S.barras}         # (22)-(23)
    th[S.slack].lb = th[S.slack].ub = 0.0                                                  # (24)

    fac = {}
    for k, L in enumerate(S.ac):
        xk = xac[k]
        bk = S.Sbase / L.x
        if isinstance(xk, (int, float)):
            if xk < 0.5:
                continue
            fac[k] = m.addVar(-L.fmax, L.fmax, name=f"fac_{tag}_{k}")                      # (12)-(13)
            m.addConstr(fac[k] == bk * (th[L.o] - th[L.r]), name=f"flujo_{tag}_{k}")       # (11)
        else:
            fac[k] = m.addVar(-L.fmax, L.fmax, name=f"fac_{tag}_{k}")
            m.addConstr(fac[k] <= L.fmax * xk, name=f"fmax_{tag}_{k}")                     # (59)
            m.addConstr(fac[k] >= -L.fmax * xk, name=f"fmin_{tag}_{k}")
            m.addConstr((1.0 / bk) * fac[k] - (th[L.o] - th[L.r]) <= (1 - xk) * M)         # (60)
            m.addConstr((1.0 / bk) * fac[k] - (th[L.o] - th[L.r]) >= -(1 - xk) * M)
    fdc = {}
    for k, L in enumerate(S.dc):
        xk = xdc[k]
        if isinstance(xk, (int, float)):
            if xk < 0.5:
                continue
            fdc[k] = m.addVar(-L.fmax, L.fmax, name=f"fdc_{tag}_{k}")                      # (14)-(16), x=1
        else:
            # f = x f^a con |f^a| <= f_max  <=>  |f| <= x f_max  (region factible identica a (61)-(63))
            fdc[k] = m.addVar(-L.fmax, L.fmax, name=f"fdc_{tag}_{k}")
            m.addConstr(fdc[k] <= L.fmax * xk, name=f"dcmax_{tag}_{k}")
            m.addConstr(fdc[k] >= -L.fmax * xk, name=f"dcmin_{tag}_{k}")

    bal = {}
    for s in S.barras:
        e = gp.LinExpr()
        for i in gb[s]:
            for b in range(len(S.term[i].bloques)):
                e += g[i, b]
        for y in eb[s]:
            e += gw[y]
        for j in db[s]:
            e += r[j]
            for h in range(len(S.dems[j].bloques)):
                e -= d[j, h]
        for k, v in fac.items():
            L = S.ac[k]
            if L.o == s:
                e -= v
            if L.r == s:
                e += v
        for k, v in fdc.items():
            L = S.dc[k]
            if L.o == s:
                e -= v
            if L.r == s:
                e += v
        bal[s] = m.addConstr(e == 0.0, name=f"bal_{tag}_{s}")                              # (10)

    dmin = {}
    for j in range(len(S.dems)):
        dmin[j] = m.addConstr(gp.quicksum(d[j, h] for h in range(len(S.dems[j].bloques)))
                              >= S.dmin_frac * _dmax_tot(S, cs, j), name=f"dmin_{tag}_{j}")   # (21)

    sw = gp.LinExpr()
    for (j, h), v in d.items():
        sw += S.dems[j].bloques[h][1] * v
    for (i, b), v in g.items():
        sw -= S.term[i].bloques[b][1] * v
    rac = gp.LinExpr()
    for j, v in r.items():
        rac += S.dems[j].cU * v
    sw -= rac                                                                              # (9)
    return dict(g=g, gw=gw, d=d, r=r, th=th, fac=fac, fdc=fdc, bal=bal, dmin=dmin, sw=sw, rac=rac)


# ----------------------------------------------------------------------
# Nivel inferior: dual (29)-(54) con (64)-(70)
# ----------------------------------------------------------------------
def add_dual_h(m, S: SistemaH, cs: Caso, xac: dict, xdc: dict, tag: str, gamma="paper"):
    GA = S.Gamma if gamma == "paper" else gamma
    GD = S.GammaDC if gamma == "paper" else gamma
    gb, eb, db = _idx_bus(S)
    INF = GRB.INFINITY

    lam = {s: m.addVar(-INF, INF, name=f"lam_{tag}_{s}") for s in S.barras}                # (39)
    vphi = {(i, b): m.addVar(0.0, INF, name=f"vphi_{tag}_{i}_{b}")
            for i, G in enumerate(S.term) for b in range(len(G.bloques))}                   # (46)
    vphiw = {y: m.addVar(0.0, INF, name=f"vphiw_{tag}_{y}") for y in range(len(S.eol))}     # (47)
    beta = {(j, h): m.addVar(0.0, INF, name=f"beta_{tag}_{j}_{h}")
            for j, D in enumerate(S.dems) for h in range(len(D.bloques))}                   # (48)
    alpha = {j: m.addVar(0.0, INF, name=f"alpha_{tag}_{j}") for j in range(len(S.dems))}    # (49)
    rho = {j: m.addVar(-INF, 0.0, name=f"rho_{tag}_{j}") for j in range(len(S.dems))}       # (50)
    ximax = {s: m.addVar(0.0, INF, name=f"ximax_{tag}_{s}") for s in S.barras if s != S.slack}   # (51)
    ximin = {s: m.addVar(-INF, 0.0, name=f"ximin_{tag}_{s}") for s in S.barras if s != S.slack}  # (52)
    chi = m.addVar(-INF, INF, name=f"chi_{tag}")                                            # (53)

    phi, pmax, pmin, pbar = {}, {}, {}, {}
    for k, L in enumerate(S.ac):
        xk = xac[k]
        if isinstance(xk, (int, float)) and xk < 0.5:
            continue
        if isinstance(xk, (int, float)):
            lo, hi = (-GA, GA) if GA is not None else (-INF, INF)
            phi[k] = m.addVar(lo, hi, name=f"phi_{tag}_{k}")
        else:
            phi[k] = m.addVar(-INF, INF, name=f"phi_{tag}_{k}")                             # (40)
        pmax[k] = m.addVar(0.0, INF, name=f"pmax_{tag}_{k}")                                # (41)
        pmin[k] = m.addVar(-INF, 0.0, name=f"pmin_{tag}_{k}")                               # (42)
    pdc, pdmax, pdmin, pdbar = {}, {}, {}, {}
    for k, L in enumerate(S.dc):
        xk = xdc[k]
        if isinstance(xk, (int, float)) and xk < 0.5:
            continue
        if isinstance(xk, (int, float)):
            lo, hi = (-GD, GD) if GD is not None else (-INF, INF)
            pdc[k] = m.addVar(lo, hi, name=f"pdc_{tag}_{k}")
        else:
            pdc[k] = m.addVar(-INF, INF, name=f"pdc_{tag}_{k}")                             # (43)
        pdmax[k] = m.addVar(0.0, INF, name=f"pdmax_{tag}_{k}")                              # (44)
        pdmin[k] = m.addVar(-INF, 0.0, name=f"pdmin_{tag}_{k}")                             # (45)

    for i, G in enumerate(S.term):                                                          # (30)
        for b, (_, lg) in enumerate(G.bloques):
            m.addConstr(lam[G.bus] + vphi[i, b] >= -lg)
    for y, E in enumerate(S.eol):                                                           # (31)
        m.addConstr(lam[E.bus] + vphiw[y] >= 0.0)
    for j, D in enumerate(S.dems):
        for h, (_, ld) in enumerate(D.bloques):                                             # (32)
            m.addConstr(-lam[D.bus] + beta[j, h] + rho[j] >= ld)
        m.addConstr(lam[D.bus] + alpha[j] >= -D.cU)                                         # (33)
    for k in phi:                                                                           # (34)
        L = S.ac[k]
        m.addConstr(-lam[L.o] + lam[L.r] + phi[k] + pmax[k] + pmin[k] == 0.0)

    xphi = {}                                                                               # x_k * phi_k
    for k in phi:
        xk = xac[k]
        if isinstance(xk, (int, float)):
            xphi[k] = phi[k]
        else:
            pbar[k] = m.addVar(-INF, INF, name=f"pbar_{tag}_{k}")                           # (69)-(70)
            xphi[k] = phi[k] - pbar[k]
            m.addConstr(phi[k] - pbar[k] <= xk * GA)
            m.addConstr(phi[k] - pbar[k] >= -xk * GA)
            m.addConstr(pbar[k] <= (1 - xk) * GA)
            m.addConstr(pbar[k] >= -(1 - xk) * GA)
    for s in S.barras:                                                                      # (37)-(38) -> (67)-(68)
        e = gp.LinExpr()
        for k in phi:
            L = S.ac[k]
            bk = S.Sbase / L.x
            if L.o == s:
                e -= bk * xphi[k]
            if L.r == s:
                e += bk * xphi[k]
        e += chi if s == S.slack else (ximax[s] + ximin[s])
        m.addConstr(e == 0.0)

    for k in pdc:                                                                           # (35)-(36) -> (64)-(66)
        L = S.dc[k]
        xk = xdc[k]
        m.addConstr(-lam[L.o] + lam[L.r] + pdc[k] == 0.0)                                   # (35)
        if isinstance(xk, (int, float)):
            xpd = pdc[k]
        else:
            pdbar[k] = m.addVar(-INF, INF, name=f"pdbar_{tag}_{k}")
            xpd = pdc[k] - pdbar[k]
            m.addConstr(pdc[k] - pdbar[k] <= xk * GD)                                       # (65)
            m.addConstr(pdc[k] - pdbar[k] >= -xk * GD)
            m.addConstr(pdbar[k] <= (1 - xk) * GD)                                          # (66)
            m.addConstr(pdbar[k] >= -(1 - xk) * GD)
        m.addConstr(-xpd + pdmax[k] + pdmin[k] == 0.0)                                      # (36)/(64)

    obj = gp.LinExpr()                                                                      # (29)
    for k in phi:
        obj += (pmax[k] - pmin[k]) * S.ac[k].fmax
    for k in pdc:
        obj += (pdmax[k] - pdmin[k]) * S.dc[k].fmax
    for (i, b), v in vphi.items():
        obj += v * S.term[i].bloques[b][0]
    for y, v in vphiw.items():
        obj += v * _pw(S, cs, y)
    for (j, h), v in beta.items():
        obj += v * _dmax_blq(S, cs, j, h)
    for j in range(len(S.dems)):
        obj += alpha[j] * _dmax_tot(S, cs, j) + rho[j] * S.dmin_frac * _dmax_tot(S, cs, j)
    for s in ximax:
        obj += math.pi * (ximax[s] - ximin[s])
    return dict(lam=lam, phi=phi, pdc=pdc, obj=obj)


# ----------------------------------------------------------------------
# Resultado
# ----------------------------------------------------------------------
@dataclass
class ResultadoH:
    estado: str = ""
    obj: float = float("nan")                  # M$/ano (ec. (1))
    plan_ac: Dict[Tuple[int, int], int] = field(default_factory=dict)
    plan_dc: Dict[Tuple[int, int], int] = field(default_factory=dict)
    inv_total: float = float("nan")            # M$
    inv_ac: float = float("nan")
    inv_dc: float = float("nan")
    esw: float = float("nan")                  # esperanza del bienestar social, M$/ano
    cvar: float = float("nan")                 # M$/ano (solo beta > 0)
    sw_casos: List[float] = field(default_factory=list)   # $/h por caso
    t: float = float("nan")
    gap: float = float("nan")
    nodos: float = float("nan")
    max_phi: float = float("nan")

    def plan_txt(self):
        return plan_txt(self.plan_ac) + ("  |  HVDC: " + plan_txt(self.plan_dc) if self.plan_dc else "")



# ----------------------------------------------------------------------
# Solver alternativo: HiGHS (libre). Se activa con la variable de entorno
# TNEP_SOLVER=highs. El modelo se construye igual con gurobipy (sin resolverlo),
# se exporta a MPS y lo resuelve HiGHS; no requiere licencia de Gurobi completa.
# ----------------------------------------------------------------------
def _resolver_highs(m, time_limit=None, mip_gap=0.0):
    import os
    import tempfile
    import highspy
    m.update()
    fd, ruta = tempfile.mkstemp(suffix=".mps")
    os.close(fd)
    try:
        m.write(ruta)
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("mip_rel_gap", float(mip_gap))
        h.setOptionValue("mip_abs_gap", 1e-9)
        h.setOptionValue("mip_feasibility_tolerance", 1e-9)
        h.setOptionValue("primal_feasibility_tolerance", 1e-8)
        h.setOptionValue("dual_feasibility_tolerance", 1e-8)
        if time_limit:
            h.setOptionValue("time_limit", float(time_limit))
        h.readModel(ruta)
        h.run()
    finally:
        os.remove(ruta)
    st = h.getModelStatus()
    estado = {highspy.HighsModelStatus.kOptimal: "OPTIMO",
              highspy.HighsModelStatus.kInfeasible: "INFACTIBLE",
              highspy.HighsModelStatus.kTimeLimit: "LIMITE_TIEMPO"}.get(st, str(st))
    info = h.getInfo()
    vals = None
    if info.primal_solution_status == 2:            # solucion factible disponible
        x = list(h.getSolution().col_value)
        gvars = m.getVars()
        if len(x) != len(gvars):
            raise RuntimeError("HiGHS: numero de columnas distinto al del modelo")
        vals = {v.index: x[i] for i, v in enumerate(gvars)}
    return dict(estado=estado, obj=info.objective_function_value, gap=info.mip_gap, nodos=info.mip_node_count, vals=vals)


def _valor(expr, vals):
    """Valor de un Var/LinExpr de gurobipy con los valores de HiGHS."""
    if isinstance(expr, (int, float)):
        return float(expr)
    if isinstance(expr, gp.Var):
        return vals[expr.index]
    return expr.getConstant() + sum(expr.getCoeff(i) * vals[expr.getVar(i).index] for i in range(expr.size()))


# ----------------------------------------------------------------------
# MILP bilevel (55)-(58) y su equivalente de un nivel (oraculo)
# ----------------------------------------------------------------------
def resolver_h(S: SistemaH, casos: List[Caso], grupos: Optional[List[Grupo]] = None,
               alpha=0.8, beta=0.0, metodo="oraculo", fix=None, simetria=True,
               presupuesto=True, time_limit=None, mip_gap=0.0,
               M=None, gamma="paper", threads=None, verbose=False, extra=None) -> ResultadoH:
    """max (1) s.a. (2)-(8), mercado de cada caso (primal; y dual + dualidad fuerte
    si metodo='dualidad').  metodo='oraculo': MILP de un solo nivel (solo primal),
    valido porque (1) es no decreciente en cada SW_c (ver el reporte, Sec. IV-D).
    fix = (plan_ac, plan_dc) fija el plan."""
    if beta > 0 and not grupos:
        raise ValueError("beta > 0 requiere los grupos de riesgo")
    m = _nuevo_modelo("heydari", verbose, time_limit, mip_gap, threads, extra)
    fixa = fixd = None
    if fix is not None:
        fixa = _x_fijo(S.ac, fix[0], True)
        fixd = _x_fijo(S.dc, fix[1], False)
    xac, xdc = {}, {}
    for k, L in enumerate(S.ac):
        if L.existente:
            xac[k] = 1.0
        else:
            v = m.addVar(0.0, 1.0, vtype=GRB.BINARY, name=f"xac_{k}")
            if fixa is not None:
                v.lb = v.ub = fixa[k]
            xac[k] = v
    for k, L in enumerate(S.dc):
        v = m.addVar(0.0, 1.0, vtype=GRB.BINARY, name=f"xdc_{k}")
        if fixd is not None:
            v.lb = v.ub = fixd[k]
        xdc[k] = v
    if simetria and fix is None:
        for lineas, x, cand in ((S.ac, xac, S.cand_ac), (S.dc, xdc, S.cand_dc)):
            por = {}
            for k in cand:
                por.setdefault((lineas[k].o, lineas[k].r), []).append(k)
            for ks in por.values():
                for a, b in zip(ks[:-1], ks[1:]):
                    m.addConstr(x[a] >= x[b])
    inv = (gp.quicksum(S.ac[k].costo * xac[k] for k in S.cand_ac)
           + gp.quicksum(S.dc[k].costo * xdc[k] for k in S.cand_dc))
    if presupuesto:
        m.addConstr(inv <= S.cmax, name="presupuesto")                                       # (3)

    prim, dual = [], []
    for c, cs in enumerate(casos):
        P = add_primal_h(m, S, cs, xac, xdc, f"c{c}", M=M)
        prim.append(P)
        if metodo == "dualidad":
            D = add_dual_h(m, S, cs, xac, xdc, f"c{c}", gamma=gamma)
            dual.append(D)
            m.addConstr(P["sw"] == D["obj"], name=f"dual_fuerte_{c}")                        # (54)

    inv_anual = S.crf * inv
    esw = gp.quicksum((cs.horas / 1e6) * prim[c]["sw"] for c, cs in enumerate(casos))
    obj = esw - inv_anual
    cvar_expr = None
    if beta > 0:
        eta = m.addVar(-GRB.INFINITY, GRB.INFINITY, name="eta")
        Sw = [m.addVar(0.0, GRB.INFINITY, name=f"S_{w}") for w in range(len(grupos))]       # (8)
        for w, gr in enumerate(grupos):
            beneficio = gp.quicksum((hrs / 1e6) * prim[c]["sw"] for c, hrs in gr.miembros) - inv_anual
            m.addConstr(eta - beneficio <= Sw[w], name=f"riesgo_{w}")                        # (2)
        cvar_expr = eta - (1.0 / (1.0 - alpha)) * gp.quicksum(gr.delta * Sw[w] for w, gr in enumerate(grupos))
        obj = obj + beta * cvar_expr                                                         # (1)
    m.setObjective(obj, GRB.MAXIMIZE)

    t0 = time.time()
    import os
    if os.environ.get("TNEP_SOLVER", "gurobi").lower() == "highs":
        sol = _resolver_highs(m, time_limit, mip_gap)
        res = ResultadoH(estado=sol["estado"], t=time.time() - t0)
        V = sol["vals"]
        if V is not None:
            xa = {k: (v if isinstance(v, float) else V[v.index]) for k, v in xac.items()}
            xd = {k: V[v.index] for k, v in xdc.items()}
            res.plan_ac, res.plan_dc = _a_plan(S.ac, xa, True), _a_plan(S.dc, xd, False)
            res.inv_total, res.inv_ac, res.inv_dc = costo_plan(S, res.plan_ac, res.plan_dc)
            res.obj = sol["obj"]
            res.sw_casos = [_valor(prim[c]["sw"], V) for c in range(len(casos))]
            res.esw = sum(cs.horas * res.sw_casos[c] for c, cs in enumerate(casos)) / 1e6
            if cvar_expr is not None:
                res.cvar = _valor(cvar_expr, V)
            res.gap, res.nodos = sol["gap"], sol["nodos"]
        return res
    m.optimize()
    res = ResultadoH(estado=_estado(m), t=time.time() - t0)
    if m.SolCount > 0:
        xa = {k: (v.X if not isinstance(v, float) else v) for k, v in xac.items()}
        xd = {k: v.X for k, v in xdc.items()}
        res.plan_ac, res.plan_dc = _a_plan(S.ac, xa, True), _a_plan(S.dc, xd, False)
        res.inv_total, res.inv_ac, res.inv_dc = costo_plan(S, res.plan_ac, res.plan_dc)
        res.obj = m.ObjVal
        res.sw_casos = [prim[c]["sw"].getValue() for c in range(len(casos))]
        res.esw = sum(cs.horas * res.sw_casos[c] for c, cs in enumerate(casos)) / 1e6
        if cvar_expr is not None:
            res.cvar = cvar_expr.getValue()
        if dual:
            res.max_phi = max([abs(v.X) for D in dual for v in D["phi"].values()]
                              + [abs(v.X) for D in dual for v in D["pdc"].values()] + [0.0])
        try:
            res.gap, res.nodos = m.MIPGap, m.NodeCount
        except gp.GurobiError:
            pass
    return res


# ----------------------------------------------------------------------
# Evaluacion de un plan fijo (LP de mercado de cada caso)
# ----------------------------------------------------------------------
def cvar_discreto(beneficios: List[float], deltas: List[float], alpha: float) -> float:
    """CVaR_alpha de beneficios con probabilidades delta: max_eta eta - 1/(1-alpha) E[(eta-p)+]
    (Rockafellar-Uryasev; el maximo se alcanza en un punto de la muestra)."""
    best = -float("inf")
    for eta in beneficios:
        v = eta - sum(d * max(eta - p, 0.0) for d, p in zip(deltas, beneficios)) / (1.0 - alpha)
        best = max(best, v)
    return best


def evaluar_h(S: SistemaH, plan_ac: dict, plan_dc: dict, casos: List[Caso],
              grupos: Optional[List[Grupo]] = None, alpha=0.8, beta=0.0, detalle=False):
    """Resuelve el LP de mercado de cada caso con el plan fijo. Devuelve un dict con
    esw (M$/ano), obj (ec. (1)), cvar, sw_casos y, con detalle=True, energias (GWh/ano)
    y estadisticas de LMP."""
    xa = _x_fijo(S.ac, plan_ac, True)
    xd = _x_fijo(S.dc, plan_dc, False)
    sws, det = [], []
    for c, cs in enumerate(casos):
        m = _nuevo_modelo("lp_h")
        P = add_primal_h(m, S, cs, xa, xd, "p")
        m.setObjective(P["sw"], GRB.MAXIMIZE)
        m.optimize()
        if m.Status != GRB.OPTIMAL:
            raise RuntimeError(f"LP de mercado no optimo (caso {cs.nombre}): {m.Status}")
        sws.append(m.ObjVal)
        if detalle:
            det.append(dict(g=sum(v.X for v in P["g"].values()), gw=sum(v.X for v in P["gw"].values()),
                            d=sum(v.X for v in P["d"].values()), r=sum(v.X for v in P["r"].values()),
                            lmp=[-P["bal"][s].Pi for s in S.barras]))
    tot, ca, cd = costo_plan(S, plan_ac, plan_dc)
    esw = sum(cs.horas * sws[c] for c, cs in enumerate(casos)) / 1e6
    out = dict(esw=esw, inv_total=tot, inv_ac=ca, inv_dc=cd, sw_casos=sws, cvar=float("nan"))
    obj = esw - S.crf * tot
    if beta > 0 and grupos:
        ben = [sum(hrs * sws[c] for c, hrs in gr.miembros) / 1e6 - S.crf * tot for gr in grupos]
        out["cvar"] = cvar_discreto(ben, [gr.delta for gr in grupos], alpha)
        obj += beta * out["cvar"]
        out["beneficios"] = ben
    out["obj"] = obj
    if detalle:
        h = [cs.horas for cs in casos]
        out["GWh_termica"] = sum(h[c] * det[c]["g"] for c in range(len(casos))) / 1e3
        out["GWh_eolica"] = sum(h[c] * det[c]["gw"] for c in range(len(casos))) / 1e3
        out["GWh_consumo"] = sum(h[c] * det[c]["d"] for c in range(len(casos))) / 1e3     # incluye lo racionado
        out["GWh_racion"] = sum(h[c] * det[c]["r"] for c in range(len(casos))) / 1e3
        out["det"] = det
        # LMP: el paper no define sus estadisticos; la definicion que reproduce las filas
        # "Average/Standard deviation of LMPs" de su Tabla 4 (casos 1-3, 4-6 cifras) es la
        # media y desviacion tipica poblacional de TODOS los pares (caso, barra), sin
        # ponderar por horas ni probabilidades.
        vals = [p for c in range(len(casos)) for p in det[c]["lmp"]]
        mu = sum(vals) / len(vals)
        out["lmp_medio"] = mu
        out["lmp_std"] = math.sqrt(sum((p - mu) ** 2 for p in vals) / len(vals))
        # variante ponderada por horas (no coincide con el paper)
        wts = [cs.horas for cs in casos for _ in S.barras]
        W = sum(wts)
        muh = sum(w * p for w, p in zip(wts, vals)) / W
        out["lmp_medio_h"] = muh
        out["lmp_std_h"] = math.sqrt(sum(w * (p - muh) ** 2 for w, p in zip(wts, vals)) / W)
    return out


def rango_energia_h(S: SistemaH, plan_ac: dict, plan_dc: dict, casos: List[Caso], tol=1e-5):
    """Rango (min, max) en GWh/ano de la generacion termica, eolica, el consumo
    (incluye lo racionado) y el racionamiento sobre TODOS los despachos igualmente
    optimos de cada caso (cuando hay empates de precio el despacho no es unico, y
    las energias publicadas dependen del desempate del solver).  Como los casos
    son independientes, el rango anual es la suma de los rangos por caso."""
    xa = _x_fijo(S.ac, plan_ac, True)
    xd = _x_fijo(S.dc, plan_dc, False)
    claves = ("termica", "eolica", "consumo", "racion")
    lo = {k: 0.0 for k in claves}
    hi = {k: 0.0 for k in claves}
    for cs in casos:
        m = _nuevo_modelo("lp_rango")
        P = add_primal_h(m, S, cs, xa, xd, "p")
        m.setObjective(P["sw"], GRB.MAXIMIZE)
        m.optimize()
        sw_opt = m.ObjVal
        m.addConstr(P["sw"] >= sw_opt - tol)
        exprs = dict(termica=gp.quicksum(P["g"].values()), eolica=gp.quicksum(P["gw"].values()),
                     consumo=gp.quicksum(P["d"].values()), racion=gp.quicksum(P["r"].values()))
        for k, e in exprs.items():
            for sentido, dest in ((GRB.MINIMIZE, lo), (GRB.MAXIMIZE, hi)):
                m.setObjective(e, sentido)
                m.optimize()
                dest[k] += cs.horas * m.ObjVal / 1e3
    return {k: (lo[k], hi[k]) for k in claves}


def dual_lp_h(S: SistemaH, plan_ac: dict, plan_dc: dict, cs: Caso, gamma=None):
    """LP dual con el plan fijo (gamma=None: sin cota). Devuelve (objetivo, max|phi|)."""
    xa = _x_fijo(S.ac, plan_ac, True)
    xd = _x_fijo(S.dc, plan_dc, False)
    m = _nuevo_modelo("lp_dual_h")
    D = add_dual_h(m, S, cs, xa, xd, "d", gamma=gamma)
    m.setObjective(D["obj"], GRB.MINIMIZE)
    m.optimize()
    if m.Status != GRB.OPTIMAL:
        return None, float("nan")
    mp = max([abs(v.X) for v in D["phi"].values()] + [abs(v.X) for v in D["pdc"].values()] + [0.0])
    return m.ObjVal, mp
