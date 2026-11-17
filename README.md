# Proyecto 4: TNEP bilevel en ambiente de mercado e incertidumbre eólica con HVDC

ILE4100 Planeamiento de Sistemas de Potencia, Universidad de los Andes.
Autores: Jaime Rueda, Fabio Silva y Andrés Méndez.

Implementación en Python + Gurobi de:

* **Reto 1.** El modelo bilevel de planeamiento de la expansión de la transmisión de Garcés, Conejo, García-Bertrand y
  Romero (IEEE Trans. Power Syst. 24(3), 2009), resuelto como MILP por dualidad fuerte y como MILP de un nivel,
  validado contra sus resultados en los sistemas de **Garver** e **IEEE RTS de 24 barras** (planes, inversiones y
  bienestar medio, incluidas contingencias y la multicontingencia de 22 escenarios).
* **Reto 2.** La extensión de Heydari y Barforoushi (Electr. Power Syst. Res. 229, 110182, 2024) con líneas HVDC,
  generación eólica y CVaR, validada en el caso de Garver (Tabla 4 del artículo), con tres metodologías de
  incertidumbre: segmentación por curva de duración, Point Estimate Method de Hong (1998) y simulación Monte Carlo con
  reducción de escenarios por k-means.

Los resultados y su análisis están en el informe IEEE (`reporte/`).

## Contenido

| Carpeta | Qué es |
|---|---|
| `codigo/` | Código fuente (modelo, datos, experimentos y verificaciones) |
| `resultados/` | Reportes `.md` y `.json` que produce cada script de `codigo/` |
| `reporte/` | Fuente LaTeX del informe IEEE, con las tablas y figuras y el script que las genera |

Archivos principales de `codigo/`: `sistema.py` y `tnep.py` (modelo de Garcés), `tnep_heydari.py` (HVDC, eólica, CVaR),
`incertidumbre.py` (segmentación, PEM y Monte Carlo + k-means), `datos_garver.py`, `datos_rts24.py` y `datos_heydari.py` (datos y resultados
publicados), `reto1_*.py`, `reto2_garver.py` y `reto2_semillas.py` (experimentos) y `verificaciones_*.py`.

## Requisitos

Python 3.14 (probado) y Gurobi 13.0.3 con licencia (la académica basta; **no** se incluye ninguna licencia).

Sin licencia completa de Gurobi, los modelos de Heydari (Reto 2) se pueden resolver con el solver libre **HiGHS**
(`pip install highspy`) definiendo la variable de entorno `TNEP_SOLVER=highs`: el modelo se arma con gurobipy (sin
resolverlo), se exporta a MPS y lo resuelve HiGHS. HiGHS reproduce los óptimos de Gurobi de la Tabla 4 (12 MILP, diferencia
<= 2e-10). La metodología 3 se corrió así.

    pip install -r requirements.txt

## Reproducir

Los scripts se ejecutan desde `codigo/` y escriben en `../resultados/`.

    cd codigo

    # Reto 1
    python verificaciones_datos.py          # costos de los planes publicados y consistencia de tablas
    python verificaciones_modelo.py         # primal = dual
    python reto1_garver.py                  # Garver: Tablas IV y V
    python reto1_rts24.py sin               # RTS-24: Tabla VIII
    python reto1_rts24.py simple            # RTS-24: Cuadro 6.9 (contingencia simple)
    python reto1_rts24.py multi 7200        # RTS-24: Tabla X (22 escenarios; la dualidad completa tarda ~35 min)

    # Reto 2
    python verificaciones_heydari.py        # 36 comprobaciones
    python reto2_garver.py seg              # Tabla 4, metodologia 1 (segmentacion)
    python reto2_garver.py dual 3600        # formulacion de dualidad
    python reto2_garver.py pem              # metodologia 2 (PEM)
    python reto2_garver.py mc               # metodologia 3 (Monte Carlo + k-means; ~11 min con HiGHS)
    python reto2_semillas.py                # metodologia 3: K = 20 con 5 semillas (casos 2-4; ~15 min con HiGHS)
    #   con HiGHS en Windows:  set TNEP_SOLVER=highs   (Linux/Mac: TNEP_SOLVER=highs python ...)
    python reto2_garver.py sens             # sensibilidad al coeficiente de demanda

Tiempos de referencia: Ryzen 7 5800H, 16 hilos.

`digitalizar_fig2_heydari.py` lee la Fig. 2 del artículo de Heydari, que **no** se incluye por derechos de autor:
colocar el PDF (`1-s2.0-S0378779624000713-main.pdf`) en la raíz del repositorio o pasar su ruta como argumento.

## Datos y su procedencia

Cada dato es publicado, tomado de una fuente citada o **inferido**; las inferencias se declaran en el informe
(Sección III-B): los niveles de viento y la curva de duración del Garver de Heydari (digitalizados de la Fig. 2, sesgo
de +0.6 % en energía eólica) y los bloques exactos y pesos 1/28 del RTS-24. Las líneas HVAC/HVDC del Garver de Heydari
vienen del repositorio GP de la UTP (http://academia.utp.edu.co/planeamiento/?p=2521); una copia de su tabla está en
`codigo/fuente_utp_garver_hvac_hvdc.csv` y la comprueba `verificaciones_heydari.py`.

## Fuentes (no incluidas)

1. L. P. Garcés, A. J. Conejo, R. García-Bertrand, R. Romero, "A bilevel approach to transmission expansion planning
   within a market environment," IEEE Trans. Power Syst., vol. 24, no. 3, pp. 1513-1522, 2009.
2. L. P. Garcés Negrete, tesis doctoral, UNESP, Ilha Solteira, 2010.
3. R. Heydari, T. Barforoushi, "HVDC/HVAC transmission network expansion planning in electricity markets in the
   presence of wind resources," Electr. Power Syst. Res., vol. 229, 110182, 2024.
4. H. P. Hong, "An efficient point estimate method for probabilistic analysis," Reliab. Eng. Syst. Saf., vol. 59,
   no. 3, pp. 261-267, 1998.

## Alcance

Las tres metodologías del Reto 2. El enunciado nombra dos; como tercera se eligió Monte Carlo + k-means
(declarado en el informe). No se incluye el sistema del sur de Brasil de Heydari.

## Licencia

Trabajo académico del curso ILE4100. Todos los derechos reservados por los autores.
