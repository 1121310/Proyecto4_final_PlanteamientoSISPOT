# -*- coding: utf-8 -*-
"""
Modelado de la incertidumbre (eolica) para el TNEP bilevel de Heydari.

  1) segmentacion(...)  : curva de duracion de carga en bloques + 6 niveles de
                          viento por bloque con pesos iguales (metodo del paper).
  3) montecarlo_kmeans(...): simulacion Monte Carlo del viento (distribucion Beta
                          por bloque ajustada por momentos a los 6 niveles) y
                          reduccion de escenarios con k-means (metodologia 3,
                          eleccion propia: el enunciado no la nombra).
  2) pem_casos(...)     : Point Estimate Method de Hong (1998). Para cada bloque
                          de la curva de duracion se toma la intensidad eolica
                          como UNA variable aleatoria (m = 1) cuya distribucion
                          son los 6 niveles equiprobables del paper; el PEM la
                          sustituye por 2 puntos (esquema 2m) o 3 puntos
                          (esquema 2m+1).

Formulas del PEM (condiciones de momentos). El articulo original de Hong (1998)
NO se consulto: el esquema 2m+1 esta tomado de arXiv:1908.00710, ec. (30)-(33)
(que describe el metodo de Hong) y el de 2 puntos se deriva de las mismas
condiciones de momentos; ambos se autoverifican numericamente abajo.
Variable estandarizada xi = (x - mu)/sigma,
lambda_3 = asimetria, lambda_4 = curtosis:
  2m  (m = 1):  p1 + p2 = 1, p1 xi1 + p2 xi2 = 0, p1 xi1^2 + p2 xi2^2 = 1,
                p1 xi1^3 + p2 xi2^3 = lambda_3
                xi_{1,2} = lambda_3/2 +- sqrt(1 + lambda_3^2/4),
                p1 = -xi2/(xi1 - xi2),  p2 = xi1/(xi1 - xi2)
  2m+1 (m = 1): ademas p3 en xi3 = 0 y p1 xi1^4 + p2 xi2^4 = lambda_4
                xi_{1,2} = lambda_3/2 +- sqrt(lambda_4 - 3 lambda_3^2/4),
                p1 = 1/(xi1 (xi1 - xi2)),  p2 = -1/(xi2 (xi1 - xi2)),  p3 = 1 - p1 - p2
Cada esquema se verifica numericamente contra sus condiciones de momentos.

Como el bienestar anual es aditivo en los bloques, aplicar el esquema 2m+1
bloque a bloque (m = 1) da la misma esperanza que el esquema 2m+1 de Hong con
m = 5 variables independientes (una por bloque): al agrupar los puntos "en la
media" de los demas bloques, el peso de cada punto del bloque t es el de m = 1
y el de su punto central es 1 - p1 - p2. (Con m = 5 el peso w0 de Hong puede ser
negativo; agrupado por bloque no lo es.) El esquema de 2 puntos (2m) se aplica
por bloque con m = 1 como un PEM independiente, no es equivalente al de m = 5.
"""
import math
from typing import List, Tuple

from tnep_heydari import Caso, Grupo


# ----------------------------------------------------------------------
# 1) segmentacion por curva de duracion (metodo de Heydari)
# ----------------------------------------------------------------------
def segmentacion(n_h, carga, viento, con_viento=True):
    """Casos y grupos de riesgo del paper. con_viento=False (caso 1): un caso por bloque.
    Escenario de riesgo omega = el omega-esimo nivel de viento de cada bloque."""
    casos, grupos = [], []
    if not con_viento:
        for t, nh in enumerate(n_h):
            casos.append(Caso(f"B{t + 1}", t, float(nh), carga[t], 0.0))
        return casos, grupos
    nw = len(viento[0])
    idx = {}
    for t, nh in enumerate(n_h):
        for k in range(nw):
            idx[t, k] = len(casos)
            casos.append(Caso(f"B{t + 1}w{k + 1}", t, nh / nw, carga[t], viento[t][k]))
    for k in range(nw):
        grupos.append(Grupo(f"w{k + 1}", 1.0 / nw, [(idx[t, k], float(n_h[t])) for t in range(len(n_h))]))
    return casos, grupos


# ----------------------------------------------------------------------
# 2) Point Estimate Method de Hong (1998)
# ----------------------------------------------------------------------
def momentos(valores: List[float], pesos: List[float] = None):
    n = len(valores)
    pesos = pesos or [1.0 / n] * n
    mu = sum(p * v for p, v in zip(pesos, valores))
    var = sum(p * (v - mu) ** 2 for p, v in zip(pesos, valores))
    sd = math.sqrt(var)
    l3 = sum(p * (v - mu) ** 3 for p, v in zip(pesos, valores)) / sd ** 3
    l4 = sum(p * (v - mu) ** 4 for p, v in zip(pesos, valores)) / sd ** 4
    return mu, sd, l3, l4


def pem_estandar(l3: float, l4: float, esquema: str) -> List[Tuple[float, float]]:
    """Puntos (xi, p) del PEM para m = 1. esquema: '2m' | '2m+1'."""
    if esquema == "2m":
        s = math.sqrt(1.0 + l3 * l3 / 4.0)
        xi1, xi2 = l3 / 2.0 + s, l3 / 2.0 - s
        pts = [(xi1, -xi2 / (xi1 - xi2)), (xi2, xi1 / (xi1 - xi2))]
        mom = 3
    elif esquema == "2m+1":
        rad = l4 - 3.0 * l3 * l3 / 4.0
        if rad <= 0:
            raise ValueError("curtosis incompatible con el esquema 2m+1")
        s = math.sqrt(rad)
        xi1, xi2 = l3 / 2.0 + s, l3 / 2.0 - s
        p1 = 1.0 / (xi1 * (xi1 - xi2))
        p2 = -1.0 / (xi2 * (xi1 - xi2))
        pts = [(xi1, p1), (xi2, p2), (0.0, 1.0 - p1 - p2)]
        mom = 4
    else:
        raise ValueError(esquema)
    # autoverificacion de las condiciones de momentos
    obj = [1.0, 0.0, 1.0, l3, l4]
    for k in range(mom + 1):
        got = sum(p * xi ** k for xi, p in pts)
        if abs(got - obj[k]) > 1e-9:
            raise AssertionError(f"PEM {esquema}: momento {k} = {got} != {obj[k]}")
    if any(p < -1e-12 for _, p in pts):
        raise AssertionError(f"PEM {esquema}: peso negativo {pts}")
    return pts


def pem_puntos(valores: List[float], esquema: str = "2m+1", pesos: List[float] = None):
    """Puntos (x, p) del PEM para la distribucion discreta (valores, pesos)."""
    mu, sd, l3, l4 = momentos(valores, pesos)
    return [(mu + xi * sd, p) for xi, p in pem_estandar(l3, l4, esquema)]


def pem_casos(n_h, carga, viento, esquema="2m+1"):
    """Casos del PEM por bloque (sin grupos de riesgo: el CVaR del caso 4 no se
    define con pesos que dependen del bloque). Devuelve (casos, info)."""
    casos, info = [], []
    for t, nh in enumerate(n_h):
        pts = pem_puntos(viento[t], esquema)
        for k, (x, p) in enumerate(pts):
            if x < 0:
                raise ValueError(f"bloque {t + 1}: punto PEM con viento negativo {x:.4f}")
            casos.append(Caso(f"B{t + 1}p{k + 1}", t, nh * p, carga[t], x))
        info.append(dict(bloque=t + 1, momentos=momentos(viento[t]), puntos=pts))
    return casos, info


# ----------------------------------------------------------------------
# 3) Monte Carlo + reduccion de escenarios por k-means
# ----------------------------------------------------------------------
def beta_momentos(valores: List[float]):
    """Parametros (a, b) de una Beta en [0, 1] con la media y la varianza de los niveles."""
    mu, sd, _, _ = momentos(valores)
    k = mu * (1.0 - mu) / sd ** 2 - 1.0
    if k <= 0:
        raise ValueError("varianza incompatible con una Beta en [0, 1]")
    return mu * k, (1.0 - mu) * k


def _kmeans(X, K, rng, iters=300):
    """k-means (Lloyd) con inicializacion k-means++. Devuelve (centroides, etiquetas)."""
    import numpy as np
    N = X.shape[0]
    C = [X[rng.integers(N)]]
    for _ in range(1, K):
        d2 = np.min(((X[:, None, :] - np.array(C)[None, :, :]) ** 2).sum(-1), axis=1)
        C.append(X[rng.choice(N, p=d2 / d2.sum())])
    C = np.array(C)
    lab = None
    for _ in range(iters):
        nuevo = np.argmin(((X[:, None, :] - C[None, :, :]) ** 2).sum(-1), axis=1)
        if lab is not None and np.array_equal(nuevo, lab):
            break
        lab = nuevo
        for k in range(K):
            if np.any(lab == k):
                C[k] = X[lab == k].mean(axis=0)
    return C, lab


def montecarlo_kmeans(n_h, carga, viento, N=2000, K=6, semilla=2026, comonotono=True):
    """Metodologia 3. Se generan N escenarios ANUALES de viento: en cada bloque t la
    intensidad sigue una Beta(a_t, b_t) con la media y la varianza de los 6 niveles del
    paper. comonotono=True empareja las muestras de los bloques por rango (misma
    dependencia que los escenarios de riesgo del paper: el w-esimo nivel de cada
    bloque); comonotono=False las toma independientes. k-means reduce los N vectores
    (5 bloques) a K escenarios representativos (centroides) con probabilidad = fraccion
    de muestras del grupo. Cada escenario k da un mercado por bloque con n_t * p_k horas
    y un grupo de riesgo (para el CVaR del caso 4). Devuelve (casos, grupos, info)."""
    import numpy as np
    rng = np.random.default_rng(semilla)
    T = len(n_h)
    ab = [beta_momentos(v) for v in viento]
    X = np.column_stack([rng.beta(a, b, N) for a, b in ab])
    if comonotono:
        X = np.sort(X, axis=0)[::-1]                 # mismo rango en todos los bloques
        X = X[rng.permutation(N)]
    C, lab = _kmeans(X, K, rng)
    p = np.bincount(lab, minlength=K) / N
    orden = np.argsort(-C.mean(axis=1))              # del escenario mas ventoso al menos
    C, p = C[orden], p[orden]
    casos, grupos = [], []
    for k in range(K):
        if p[k] == 0:
            continue
        miembros = []
        for t in range(T):
            miembros.append((len(casos), float(n_h[t])))
            casos.append(Caso(f"B{t + 1}k{k + 1}", t, n_h[t] * float(p[k]), carga[t], float(C[k, t])))
        grupos.append(Grupo(f"k{k + 1}", float(p[k]), miembros))
    # autoverificacion: horas y probabilidades
    assert abs(sum(c.horas for c in casos) - sum(n_h)) < 1e-6
    assert abs(sum(g.delta for g in grupos) - 1.0) < 1e-12
    info = dict(N=N, K=K, semilla=semilla, comonotono=comonotono, beta=ab,
                media_muestral=X.mean(axis=0).tolist(),
                media_reducida=(p[:, None] * C).sum(axis=0).tolist(),
                media_niveles=[sum(v) / len(v) for v in viento],
                centroides=C.tolist(), prob=p.tolist())
    return casos, grupos, info
