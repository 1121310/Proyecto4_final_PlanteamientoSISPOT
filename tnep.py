# -*- coding: utf-8 -*-
"""
Modelo bilevel de planeamiento de la expansion de la transmision en un
ambiente de mercado:

    L. P. Garces, A. J. Conejo, R. Garcia-Bertrand, R. Romero,
    "A Bilevel Approach to Transmission Expansion Planning Within a Market
    Environment", IEEE Trans. Power Systems, 24(3), 1513-1522, 2009.
    (Ecuaciones (1)-(45) referenciadas en los comentarios.)

Contenido
---------
add_primal   : problema de mercado del nivel inferior, ec. (5)-(19), con la
               linealizacion (40)-(41) cuando x_k es variable binaria.
add_dual     : problema dual (21)-(34) con la linealizacion (42)-(45).
resolver_bilevel : MILP equivalente (36)-(39):  max (1)  s.a. (2)-(4),
               primal, dual y dualidad fuerte (35) por escenario.
               metodo="dualidad" -> formulacion del paper
               metodo="oraculo"  -> MILP de un solo nivel con solo el primal
                                    (valido porque (1) alinea ambos niveles,
                                    Sec. II-B). Sirve como verificacion.
evaluar_plan : bienestar social con el plan FIJO (LP del nivel inferior).
dual_lp      : LP dual con el plan fijo (verificacion de dualidad fuerte).
resolver_costmin : modelo de minimizacion de costos (un escenario a la vez).
"""
import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import gurobipy as gp
from gurobipy import GRB

from sistema import Sistema, Escenario, plan_a_texto

M_EUR = 1.0e6     # EUR por M EUR


# ----------------------------------------------------------------------
# utilidades
# ----------------------------------------------------------------------
_ENV = None


def _env() -> gp.Env:
    """Un unico entorno Gurobi compartido (evita re-licenciar en cada modelo)."""
    global _ENV
    if _ENV is None:
        e = gp.Env(empty=True)
        e.setParam("OutputFlag", 0)
        e.start()
        _ENV = e
    return _ENV


def _nuevo_modelo(nombre: str, verbose=False, time_limit=None, mip_gap=0.0,
                  threads=None, extra=None) -> gp.Model:
    m = gp.Model(nombre, env=_env())
    m.Params.OutputFlag = 1 if verbose else 0
    # M = 1000 en las desigualdades disyuntivas -> tolerancia entera estricta
    m.Params.IntFeasTol = 1e-9
    m.Params.FeasibilityTol = 1e-8
    m.Params.OptimalityTol = 1e-8
    m.Params.MIPGap = mip_gap
    m.Params.MIPGapAbs = 1e-9
    if time_limit:
        m.Params.TimeLimit = time_limit
    if threads:
        m.Params.Threads = threads
    for k, v in (extra or {}).items():
        m.setParam(k, v)
    return m


def plan_a_x(S: Sistema, plan: Dict[Tuple[int, int], int]) -> Dict[int, float]:
    """Plan {(o,r): n} -> x_k (0/1) para todas las lineas (existentes = 1).
    Las n primeras copias candidatas de cada corredor se consideran construidas
    (las copias de un corredor son identicas)."""
    x = {}
    contador: Dict[Tuple[int, int], int] = {}
    for k, L in enumerate(S.lineas):
        if L.existente:
            x[k] = 1.0
        else:
            c = contador.get((L.o, L.r), 0)
            x[k] = 1.0 if c < plan.get((L.o, L.r), 0) else 0.0
            contador[(L.o, L.r)] = c + 1
    for (o, r), n in plan.items():
        cand = sum(1 for L in S.lineas if (L.o, L.r) == (o, r) and not L.existente)
        if n > cand:
            raise ValueError(f"plan pide {n} circuitos en {o}-{r}; maximo {cand}")
    return x


def x_a_plan(S: Sistema, xval: Dict[int, float]) -> Dict[Tuple[int, int], int]:
    plan: Dict[Tuple[int, int], int] = {}
    for k in S.candidatas:
        if xval[k] > 0.5:
            c = (S.lineas[k].o, S.lineas[k].r)
            plan[c] = plan.get(c, 0) + 1
    return plan


def _indices_por_bus(S: Sistema):
    gen_bus = {s: [] for s in S.barras}
    dem_bus = {s: [] for s in S.barras}
    for i, g in enumerate(S.gens):
        gen_bus[g.bus].append(i)
    for j, d in enumerate(S.dems):
        dem_bus[d.bus].append(j)
    return gen_bus, dem_bus


# ----------------------------------------------------------------------
# Nivel inferior: problema primal  (5)-(19), linealizado con (40)-(41)
# ----------------------------------------------------------------------
def add_primal(m: gp.Model, S: Sistema, sc: Escenario, x: dict, tag: str, M=None):
    """x[k] = float (0/1, plan fijo: red exacta) o gp.Var (MILP, ec. (40)-(41)).
    Devuelve dict con variables, expresion de bienestar social (5) y
    expresion del costo de racionamiento."""
    M = S.M if M is None else M
    gen_bus, dem_bus = _indices_por_bus(S)
    fuera = S.indices_fuera(sc)
    disp = [k for k in range(len(S.lineas)) if k not in fuera]

    g = {(i, b): m.addVar(0.0, S.gens[i].bloques[b][0], name=f"g_{tag}_{i}_{b}")
         for i, G in enumerate(S.gens) for b in range(len(G.bloques))}                 # (10),(17)
    d = {(j, h): m.addVar(0.0, S.dem_max_bloque(sc, j, h), name=f"d_{tag}_{j}_{h}")
         for j, D in enumerate(S.dems) for h in range(len(D.bloques))}                 # (11),(19)
    r = {j: m.addVar(0.0, S.dem_max(sc, j), name=f"r_{tag}_{j}")
         for j in range(len(S.dems))}                                                 # (12),(18)
    th = {s: m.addVar(-math.pi, math.pi, name=f"th_{tag}_{s}") for s in S.barras}    # (14),(15)
    th[S.slack].lb = 0.0
    th[S.slack].ub = 0.0                                                              # (16)

    f = {}
    for k in disp:
        L = S.lineas[k]
        xk = x[k]
        bk = S.Sbase / L.x                       # MW/rad
        if isinstance(xk, (int, float)):
            if xk < 0.5:
                continue                         # linea no construida: no existe en la red
            f[k] = m.addVar(-L.fmax, L.fmax, name=f"f_{tag}_{k}")                       # (8),(9)
            m.addConstr(f[k] == bk * (th[L.o] - th[L.r]), name=f"flujo_{tag}_{k}")      # (7)
        else:
            f[k] = m.addVar(-L.fmax, L.fmax, name=f"f_{tag}_{k}")
            m.addConstr(f[k] <= L.fmax * xk, name=f"fmax_{tag}_{k}")                    # (40)
            m.addConstr(f[k] >= -L.fmax * xk, name=f"fmin_{tag}_{k}")
            # (41):  -(1-x)M <= f/b - (th_o - th_r) <= (1-x)M
            m.addConstr((1.0 / bk) * f[k] - (th[L.o] - th[L.r]) <= (1 - xk) * M, name=f"dj1_{tag}_{k}")
            m.addConstr((1.0 / bk) * f[k] - (th[L.o] - th[L.r]) >= -(1 - xk) * M, name=f"dj2_{tag}_{k}")

    bal = {}
    for s in S.barras:
        e = gp.LinExpr()
        for i in gen_bus[s]:
            for b in range(len(S.gens[i].bloques)):
                e += g[i, b]
        for j in dem_bus[s]:
            e += r[j]
            for h in range(len(S.dems[j].bloques)):
                e -= d[j, h]
        for k in f:
            L = S.lineas[k]
            if L.o == s:
                e -= f[k]
            if L.r == s:
                e += f[k]
        bal[s] = m.addConstr(e == 0.0, name=f"bal_{tag}_{s}")                           # (6)

    dmin = {}
    for j in range(len(S.dems)):
        dmin[j] = m.addConstr(gp.quicksum(d[j, h] for h in range(len(S.dems[j].bloques)))
                              >= S.dem_min(sc, j), name=f"dmin_{tag}_{j}")              # (13)

    sw = gp.LinExpr()
    for (j, h), v in d.items():
        sw += S.dems[j].bloques[h][1] * v
    for (i, b), v in g.items():
        sw -= S.gens[i].bloques[b][1] * v
    rac = gp.LinExpr()
    for j, v in r.items():
        rac += S.dems[j].cU * v
    sw -= rac                                                                           # (5)
    return dict(g=g, d=d, r=r, th=th, f=f, bal=bal, dmin=dmin, sw=sw, rac=rac, disp=disp)


# ----------------------------------------------------------------------
# Nivel inferior: problema dual (21)-(34), linealizado con (42)-(45)
# ----------------------------------------------------------------------
def add_dual(m: gp.Model, S: Sistema, sc: Escenario, x: dict, tag: str,
             gamma="paper"):
    """x[k] = float (plan fijo: dual exacto del LP reducido) o gp.Var
    (ec. (42)-(45) con Gamma_max).  gamma: 'paper' -> S.Gamma, None -> sin cota
    (solo plan fijo), o un numero."""
    G_ = S.Gamma if gamma == "paper" else gamma
    gen_bus, dem_bus = _indices_por_bus(S)
    fuera = S.indices_fuera(sc)
    disp = [k for k in range(len(S.lineas)) if k not in fuera]
    INF = GRB.INFINITY

    lam = {s: m.addVar(-INF, INF, name=f"lam_{tag}_{s}") for s in S.barras}           # (28)
    vphi = {(i, b): m.addVar(0.0, INF, name=f"vphi_{tag}_{i}_{b}")
            for i, Gn in enumerate(S.gens) for b in range(len(Gn.bloques))}           # (30)
    beta = {(j, h): m.addVar(0.0, INF, name=f"beta_{tag}_{j}_{h}")
            for j, D in enumerate(S.dems) for h in range(len(D.bloques))}             # (31)
    alpha = {j: m.addVar(0.0, INF, name=f"alpha_{tag}_{j}") for j in range(len(S.dems))}   # (32)
    rho = {j: m.addVar(-INF, 0.0, name=f"rho_{tag}_{j}") for j in range(len(S.dems))}      # (32)
    xi_max = {s: m.addVar(0.0, INF, name=f"ximax_{tag}_{s}") for s in S.barras if s != S.slack}   # (33)
    xi_min = {s: m.addVar(-INF, 0.0, name=f"ximin_{tag}_{s}") for s in S.barras if s != S.slack}  # (33)
    chi = m.addVar(-INF, INF, name=f"chi_{tag}")                                       # (34)

    phi, phimax, phimin, phibar = {}, {}, {}, {}
    for k in disp:
        L = S.lineas[k]
        xk = x[k]
        if isinstance(xk, (int, float)) and xk < 0.5:
            continue
        lo, hi = (-G_, G_) if (G_ is not None) else (-INF, INF)
        phi[k] = m.addVar(lo, hi, name=f"phi_{tag}_{k}") if isinstance(xk, (int, float)) \
            else m.addVar(-INF, INF, name=f"phi_{tag}_{k}")                            # (29)
        phimax[k] = m.addVar(0.0, INF, name=f"phimax_{tag}_{k}")                       # (29)
        phimin[k] = m.addVar(-INF, 0.0, name=f"phimin_{tag}_{k}")                      # (29)

    # (22)
    for i, Gn in enumerate(S.gens):
        for b, (_, lg) in enumerate(Gn.bloques):
            m.addConstr(lam[Gn.bus] + vphi[i, b] >= -lg, name=f"d22_{tag}_{i}_{b}")
    # (23),(24)
    for j, D in enumerate(S.dems):
        for h, (_, ld) in enumerate(D.bloques):
            m.addConstr(-lam[D.bus] + beta[j, h] + rho[j] >= ld, name=f"d23_{tag}_{j}_{h}")
        m.addConstr(lam[D.bus] + alpha[j] >= -D.cU, name=f"d24_{tag}_{j}")
    # (25)
    for k in phi:
        L = S.lineas[k]
        m.addConstr(-lam[L.o] + lam[L.r] + phi[k] + phimax[k] + phimin[k] == 0.0,
                    name=f"d25_{tag}_{k}")

    # (26)/(27) -> (42)/(43) con (44)/(45)
    xphi = {}                     # expresion "x_k * phi_k"
    for k in phi:
        xk = x[k]
        if isinstance(xk, (int, float)):
            xphi[k] = phi[k]                                  # x = 1
        else:
            phibar[k] = m.addVar(-INF, INF, name=f"phibar_{tag}_{k}")
            xphi[k] = phi[k] - phibar[k]                      # (42)-(43)
            m.addConstr(phi[k] - phibar[k] <= xk * G_, name=f"d44a_{tag}_{k}")
            m.addConstr(phi[k] - phibar[k] >= -xk * G_, name=f"d44b_{tag}_{k}")
            m.addConstr(phibar[k] <= (1 - xk) * G_, name=f"d45a_{tag}_{k}")
            m.addConstr(phibar[k] >= -(1 - xk) * G_, name=f"d45b_{tag}_{k}")
    for s in S.barras:
        e = gp.LinExpr()
        for k in phi:
            L = S.lineas[k]
            bk = S.Sbase / L.x
            if L.o == s:
                e -= bk * xphi[k]
            if L.r == s:
                e += bk * xphi[k]
        if s == S.slack:
            e += chi
        else:
            e += xi_max[s] + xi_min[s]
        m.addConstr(e == 0.0, name=f"d26_{tag}_{s}")

    # funcion objetivo dual (21)
    obj = gp.LinExpr()
    for k in phi:
        obj += (phimax[k] - phimin[k]) * S.lineas[k].fmax
    for (i, b), v in vphi.items():
        obj += v * S.gens[i].bloques[b][0]
    for (j, h), v in beta.items():
        obj += v * S.dem_max_bloque(sc, j, h)
    for j in range(len(S.dems)):
        obj += alpha[j] * S.dem_max(sc, j) + rho[j] * S.dem_min(sc, j)
    for s in xi_max:
        obj += math.pi * (xi_max[s] - xi_min[s])
    return dict(lam=lam, phi=phi, phimax=phimax, phimin=phimin, phibar=phibar,
                vphi=vphi, beta=beta, alpha=alpha, rho=rho, xi_max=xi_max,
                xi_min=xi_min, chi=chi, obj=obj)


# ----------------------------------------------------------------------
# Resultado
# ----------------------------------------------------------------------
@dataclass
class Resultado:
    estado: str = ""
    obj: float = float("nan")              # M EUR/ano
    plan: Dict[Tuple[int, int], int] = field(default_factory=dict)
    inv_total: float = float("nan")        # M EUR
    sw_prom: float = float("nan")          # M EUR/ano  (sigma * sum delta*SW)
    sw_esc: List[float] = field(default_factory=list)      # EUR/h por escenario
    rac_esc: List[float] = field(default_factory=list)     # MW racionados por escenario
    t: float = float("nan")
    gap: float = float("nan")
    nodos: float = float("nan")
    max_phi: float = float("nan")
    n_bin: int = 0
    pool: list = field(default_factory=list)   # planes empatados en el optimo (solo costmin)

    def plan_txt(self):
        return plan_a_texto(self.plan)


def _estado(m: gp.Model) -> str:
    return {GRB.OPTIMAL: "OPTIMO", GRB.TIME_LIMIT: "LIMITE_TIEMPO",
            GRB.INFEASIBLE: "INFACTIBLE", GRB.INF_OR_UNBD: "INFACT_o_NOACOT",
            GRB.UNBOUNDED: "NO_ACOTADO", GRB.SUBOPTIMAL: "SUBOPTIMO"}.get(m.Status, str(m.Status))


# ----------------------------------------------------------------------
# MILP bilevel (36)-(39)
# ----------------------------------------------------------------------
def resolver_bilevel(S: Sistema, escenarios: List[Escenario], metodo="dualidad",
                     objetivo="bienestar", plan_fijo: Optional[dict] = None,
                     presupuesto=True, simetria=True, dualidad="igualdad",
                     time_limit=None, mip_gap=0.0, verbose=False, M=None, gamma="paper",
                     threads=None, extra=None) -> Resultado:
    """Resuelve max (1) [objetivo='bienestar'] o max (20) [objetivo='racionamiento'].
    metodo='dualidad': problema (36)-(39) del paper.
    metodo='oraculo' : un solo nivel (solo primal); solo para objetivo='bienestar'.
    dualidad='igualdad' impone (35) como igualdad; 'desigualdad' impone solo
    primal >= dual (la dualidad debil da el otro sentido)."""
    if metodo == "oraculo" and objetivo != "bienestar":
        raise ValueError("el oraculo de un nivel solo es valido para el objetivo (1)")
    m = _nuevo_modelo("bilevel", verbose, time_limit, mip_gap, threads, extra)
    use_var = True
    x = {}
    fix = None if plan_fijo is None else plan_a_x(S, plan_fijo)
    for k, L in enumerate(S.lineas):
        if L.existente:
            if metodo == "dualidad":
                x[k] = m.addVar(1.0, 1.0, name=f"x_{k}")        # (3): x_k = 1
            else:
                x[k] = 1.0
        else:
            v = m.addVar(0.0, 1.0, vtype=GRB.BINARY, name=f"x_{k}")        # (4)
            if fix is not None:
                v.lb = v.ub = fix[k]
            x[k] = v
    cand = S.candidatas
    if simetria and fix is None:      # circuitos identicos de un corredor: x_copia1 >= x_copia2 >= ...
        por_cor: Dict[int, List[int]] = {}
        for k in cand:
            por_cor.setdefault(S.lineas[k].corredor, []).append(k)
        for ks in por_cor.values():
            for a, b in zip(ks[:-1], ks[1:]):
                m.addConstr(x[a] >= x[b], name=f"sim_{a}_{b}")
    if presupuesto:                                                          # (2)
        m.addConstr(gp.quicksum(S.lineas[k].costo / 1000.0 * x[k] for k in cand) <= S.cmax,
                    name="presupuesto")

    prim, dual = [], []
    for w, sc in enumerate(escenarios):
        P = add_primal(m, S, sc, x, f"w{w}", M=M)
        prim.append(P)
        if metodo == "dualidad":
            D = add_dual(m, S, sc, x, f"w{w}", gamma=gamma)
            dual.append(D)
            if dualidad == "igualdad":
                m.addConstr(P["sw"] == D["obj"], name=f"dual_fuerte_{w}")    # (35)
            else:
                m.addConstr(P["sw"] >= D["obj"], name=f"dual_fuerte_{w}")

    inv_anual = gp.quicksum(S.crf * S.lineas[k].costo / 1000.0 * x[k] for k in cand)
    if objetivo == "bienestar":
        sw_prom = gp.quicksum(sc.peso * prim[w]["sw"] for w, sc in enumerate(escenarios))
        obj = (S.sigma / M_EUR) * sw_prom - inv_anual                        # (1)
    else:
        rac_prom = gp.quicksum(sc.peso * prim[w]["rac"] for w, sc in enumerate(escenarios))
        obj = -(S.sigma / M_EUR) * rac_prom - inv_anual                      # (20)
    m.setObjective(obj, GRB.MAXIMIZE)

    t0 = time.time()
    m.optimize()
    res = Resultado(estado=_estado(m), t=time.time() - t0,
                    n_bin=sum(1 for k in cand if (fix is None)))
    if m.SolCount > 0:
        xv = {k: (x[k].X if not isinstance(x[k], float) else x[k]) for k in x}
        res.plan = x_a_plan(S, xv)
        res.obj = m.ObjVal
        res.inv_total = S.costo_plan(res.plan)
        res.sw_esc = [prim[w]["sw"].getValue() for w in range(len(escenarios))]
        res.rac_esc = [sum(v.X for v in prim[w]["r"].values()) for w in range(len(escenarios))]
        res.sw_prom = (S.sigma / M_EUR) * sum(sc.peso * res.sw_esc[w]
                                              for w, sc in enumerate(escenarios))
        if dual:
            res.max_phi = max((abs(v.X) for D in dual for v in D["phi"].values()), default=0.0)
        try:
            res.gap = m.MIPGap
            res.nodos = m.NodeCount
        except gp.GurobiError:
            pass
    return res


# ----------------------------------------------------------------------
# Evaluacion de un plan fijo (nivel inferior) y verificacion de dualidad
# ----------------------------------------------------------------------
def evaluar_plan(S: Sistema, plan: dict, escenarios: List[Escenario], detalle=False):
    """Resuelve el LP de mercado (5)-(19) de cada escenario con el plan fijo.
    Devuelve (sw_prom [M EUR/ano], [SW_w EUR/h], [MW racionados_w])."""
    x = plan_a_x(S, plan)
    sws, racs, det = [], [], []
    for w, sc in enumerate(escenarios):
        m = _nuevo_modelo("lp_plan")
        P = add_primal(m, S, sc, x, "p")
        m.setObjective(P["sw"], GRB.MAXIMIZE)
        m.optimize()
        if m.Status != GRB.OPTIMAL:
            raise RuntimeError(f"LP de mercado no optimo (escenario {sc.nombre}): {m.Status}")
        sws.append(m.ObjVal)
        racs.append(sum(v.X for v in P["r"].values()))
        if detalle:
            det.append(dict(g={k: v.X for k, v in P["g"].items()},
                            d={k: v.X for k, v in P["d"].items()},
                            r={k: v.X for k, v in P["r"].items()},
                            f={k: v.X for k, v in P["f"].items()},
                            lam={s: c.Pi for s, c in P["bal"].items()}))
    sw_prom = (S.sigma / M_EUR) * sum(sc.peso * sws[w] for w, sc in enumerate(escenarios))
    if detalle:
        return sw_prom, sws, racs, det
    return sw_prom, sws, racs


def dual_lp(S: Sistema, plan: dict, sc: Escenario, gamma=None):
    """LP dual (21)-(34) con el plan fijo. gamma=None -> sin cota en phi.
    Devuelve (objetivo dual, max|phi|)."""
    x = plan_a_x(S, plan)
    m = _nuevo_modelo("lp_dual")
    D = add_dual(m, S, sc, x, "d", gamma=gamma)
    m.setObjective(D["obj"], GRB.MINIMIZE)
    m.optimize()
    if m.Status != GRB.OPTIMAL:
        return None, float("nan")
    return m.ObjVal, max((abs(v.X) for v in D["phi"].values()), default=0.0)


def objetivo_plan(S: Sistema, plan: dict, escenarios: List[Escenario]) -> float:
    """Valor de la funcion objetivo (1) [M EUR/ano] de un plan dado:
    sigma*sum(delta*SW) - crf*inversion. Permite comprobar si el plan publicado
    es un optimo alternativo (mismo objetivo) del plan hallado."""
    sw_prom, _, _ = evaluar_plan(S, plan, escenarios)
    return sw_prom - S.crf * S.costo_plan(plan)


def resolver_centralizado(S: Sistema, sc: Escenario, frac_demanda=1.0,
                          mip_gap=0.0, verbose=False) -> Resultado:
    """Modelo 'basico y centralizado' (Garver clasico, citado en la tesis):
    min inversion  s.a. red DC, limites de generacion, demanda FIJA = frac x oferta
    total del escenario y sin racionamiento (r = 0)."""
    m = _nuevo_modelo("central", verbose, mip_gap=mip_gap)
    x = {k: (1.0 if L.existente else m.addVar(vtype=GRB.BINARY, name=f"x_{k}"))
         for k, L in enumerate(S.lineas)}
    por_cor: Dict[int, List[int]] = {}
    for k in S.candidatas:
        por_cor.setdefault(S.lineas[k].corredor, []).append(k)
    for ks in por_cor.values():
        for a, b in zip(ks[:-1], ks[1:]):
            m.addConstr(x[a] >= x[b])
    P = add_primal(m, S, sc, x, "c")
    for j in range(len(S.dems)):
        m.addConstr(gp.quicksum(P["d"][j, h] for h in range(len(S.dems[j].bloques)))
                    == frac_demanda * S.dem_max(sc, j))
        P["r"][j].ub = 0.0
    m.setObjective(gp.quicksum(S.lineas[k].costo / 1000.0 * x[k] for k in S.candidatas), GRB.MINIMIZE)
    t0 = time.time()
    m.optimize()
    res = Resultado(estado=_estado(m), t=time.time() - t0)
    if m.SolCount > 0:
        xv = {k: (x[k].X if not isinstance(x[k], float) else x[k]) for k in x}
        res.plan = x_a_plan(S, xv)
        res.obj = m.ObjVal
        res.inv_total = S.costo_plan(res.plan)
    return res


# ----------------------------------------------------------------------
# Modelo de minimizacion de costos (un escenario a la vez) - Sec. III-B
# ----------------------------------------------------------------------
def resolver_costmin(S: Sistema, sc: Escenario, ponderacion="anual", presupuesto=False,
                     simetria=True, time_limit=None, mip_gap=0.0, M=None, verbose=False,
                     threads=None, pool_n=0) -> Resultado:
    """min  inversion + costo de racionamiento  s.a. (3)-(4), (6)-(19), un escenario.
    ponderacion='anual' : crf*sum(c x)  +  sigma*sum(cU r)           [M EUR/ano]
    ponderacion='plano' : sum(c x)      +  sum(cU r)  (sin anualizar) [M EUR]"""
    m = _nuevo_modelo("costmin", verbose, time_limit, mip_gap, threads)
    x = {}
    for k, L in enumerate(S.lineas):
        x[k] = 1.0 if L.existente else m.addVar(0.0, 1.0, vtype=GRB.BINARY, name=f"x_{k}")
    cand = S.candidatas
    if simetria:
        por_cor: Dict[int, List[int]] = {}
        for k in cand:
            por_cor.setdefault(S.lineas[k].corredor, []).append(k)
        for ks in por_cor.values():
            for a, b in zip(ks[:-1], ks[1:]):
                m.addConstr(x[a] >= x[b])
    if presupuesto:
        m.addConstr(gp.quicksum(S.lineas[k].costo / 1000.0 * x[k] for k in cand) <= S.cmax)
    P = add_primal(m, S, sc, x, "c", M=M)
    inv = gp.quicksum(S.lineas[k].costo / 1000.0 * x[k] for k in cand)
    if ponderacion == "anual":
        obj = S.crf * inv + (S.sigma / M_EUR) * P["rac"]
    else:
        obj = inv + P["rac"] / M_EUR
    m.setObjective(obj, GRB.MINIMIZE)
    if pool_n:                      # enumerar optimos alternativos (empates de costo)
        m.Params.PoolSearchMode = 2
        m.Params.PoolSolutions = pool_n
        m.Params.PoolGap = 1e-9
    t0 = time.time()
    m.optimize()
    res = Resultado(estado=_estado(m), t=time.time() - t0)
    if m.SolCount > 0:
        xv = {k: (x[k].X if not isinstance(x[k], float) else x[k]) for k in x}
        res.plan = x_a_plan(S, xv)
        res.obj = m.ObjVal
        res.inv_total = S.costo_plan(res.plan)
        res.rac_esc = [sum(v.X for v in P["r"].values())]
        res.sw_esc = [P["sw"].getValue()]       # SW en el despacho (degenerado) del propio modelo
        if pool_n:
            vistos = set()
            for i in range(m.SolCount):
                m.Params.SolutionNumber = i
                if m.PoolObjVal > m.ObjVal + 1e-6:
                    continue
                xi = {k: (x[k].Xn if not isinstance(x[k], float) else x[k]) for k in x}
                pl = x_a_plan(S, xi)
                clave = tuple(sorted(pl.items()))
                if clave not in vistos:
                    vistos.add(clave)
                    res.pool.append(pl)
    return res
