# -*- coding: utf-8 -*-
"""
Estructuras de datos comunes para el modelo bilevel de expansion de la
transmision (Garces, Conejo, Garcia-Bertrand y Romero, IEEE TPWRS 24(3), 2009).

Unidades:
  - potencias en MW, precios en EUR/MWh, costos de inversion en miles de EUR
    (k EUR, tal como en las Tablas I y A.3 del paper y de la tesis).
  - reactancias en p.u. sobre una base de 100 MVA (hipotesis; el paper no la
    indica explicitamente).
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class Linea:
    o: int                # barra de envio
    r: int                # barra de recepcion
    x: float              # reactancia (p.u.)
    fmax: float           # capacidad (MW)
    costo: float          # costo de inversion de UN circuito (k EUR)
    existente: bool       # True -> ya construida (x_k = 1)
    corredor: int         # indice del corredor
    copia: int            # numero de copia dentro del corredor (0,1,2)


@dataclass
class Generador:
    bus: int
    bloques: List[Tuple[float, float]]      # (tamano MW, precio EUR/MWh)


@dataclass
class Demanda:
    bus: int
    bloques: List[Tuple[float, float]]      # (tamano MW, precio EUR/MWh)
    # costo de racionamiento = 100 x precio del primer bloque (paper, Sec. III-A)
    @property
    def cU(self) -> float:
        return 100.0 * self.bloques[0][1]


@dataclass
class Escenario:
    nombre: str
    peso: float
    factores: Dict[int, float] = field(default_factory=dict)   # bus -> factor (def. 1)
    linea_fuera: Optional[Tuple[int, int]] = None              # corredor (o,r) con contingencia

    def factor(self, bus: int) -> float:
        return self.factores.get(bus, 1.0)


@dataclass
class Sistema:
    nombre: str
    barras: List[int]
    slack: int
    lineas: List[Linea]
    gens: List[Generador]
    dems: List[Demanda]
    cmax: float                    # presupuesto de inversion (M EUR, costo TOTAL)
    sigma: float = 8760.0          # horas/ano
    crf: float = 0.10              # inversion anualizada ~ 10 % del total (paper)
    Sbase: float = 100.0           # MVA
    dmin_frac: float = 0.90        # consumo minimo = 90 % de la oferta total
    M: float = 1000.0              # constante disyuntiva (paper: 1000)
    Gamma: float = 1000.0          # cota de phi (paper: 1000)
    max_lineas_corredor: int = 3

    # ---------------- utilidades ----------------
    @property
    def candidatas(self) -> List[int]:
        return [k for k, L in enumerate(self.lineas) if not L.existente]

    @property
    def n_corredores(self) -> int:
        return max(L.corredor for L in self.lineas) + 1

    def corredor_de(self, k: int) -> Tuple[int, int]:
        L = self.lineas[k]
        return (L.o, L.r)

    def indices_fuera(self, sc: Escenario) -> set:
        """Indices de lineas NO disponibles en el escenario (contingencia).
        Sale de servicio el (primer) circuito existente del corredor; los
        circuitos candidatos del mismo corredor siguen siendo construibles."""
        if sc.linea_fuera is None:
            return set()
        for k, L in enumerate(self.lineas):
            if L.existente and (L.o, L.r) == tuple(sc.linea_fuera):
                return {k}
        raise ValueError(f"contingencia sobre corredor inexistente {sc.linea_fuera}")

    def dem_max_bloque(self, sc: Escenario, j: int, h: int) -> float:
        D = self.dems[j]
        return D.bloques[h][0] * sc.factor(D.bus)

    def dem_max(self, sc: Escenario, j: int) -> float:
        return sum(self.dem_max_bloque(sc, j, h) for h in range(len(self.dems[j].bloques)))

    def dem_min(self, sc: Escenario, j: int) -> float:
        return self.dmin_frac * self.dem_max(sc, j)

    def costo_plan(self, plan: Dict[Tuple[int, int], int]) -> float:
        """Costo TOTAL (M EUR) de un plan {(o,r): n_circuitos_nuevos}."""
        total = 0.0
        for (o, r), n in plan.items():
            cs = [L.costo for L in self.lineas
                  if not L.existente and (L.o, L.r) == (o, r)]
            if not cs:
                raise ValueError(f"corredor {(o, r)} sin candidatas")
            total += n * cs[0]
        return total / 1000.0


def construir_lineas(tabla, max_lineas=3) -> List[Linea]:
    """tabla: lista de (o, r, x, fmax, costo_kEUR, n_existentes).
    Genera n_existentes circuitos existentes y (max_lineas - n_existentes)
    candidatos por corredor."""
    lineas: List[Linea] = []
    for c, (o, r, x, fmax, costo, nex) in enumerate(tabla):
        for t in range(max_lineas):
            lineas.append(Linea(o, r, x, fmax, costo, existente=(t < nex),
                                corredor=c, copia=t))
    return lineas


def parse_plan(texto: str) -> Dict[Tuple[int, int], int]:
    """'3-5, 2-3, 4-6 (3)' -> {(3,5):1,(2,3):1,(4,6):3}"""
    plan: Dict[Tuple[int, int], int] = {}
    for tok in texto.split(","):
        tok = tok.strip()
        if not tok:
            continue
        n = 1
        if "(" in tok:
            tok, rest = tok.split("(")
            n = int(rest.strip(") "))
        a, b = tok.strip().split("-")
        plan[(int(a), int(b))] = n
    return plan


def plan_a_texto(plan: Dict[Tuple[int, int], int]) -> str:
    partes = []
    for (o, r), n in sorted(plan.items()):
        partes.append(f"{o}-{r}" + (f" ({n})" if n > 1 else ""))
    return ", ".join(partes)
