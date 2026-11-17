# -*- coding: utf-8 -*-
"""
Digitalizacion de la Fig. 2 de Heydari y Barforoushi (EPSR 229, 2024, 110182):
curva de duracion de carga (5 bloques) y 6 intensidades de viento por bloque.

El paper NO tabula estos valores (solo los grafica y remite al ref. [2],
Baringo y Conejo 2012). Aqui se leen los segmentos negros horizontales de la
imagen raster incrustada en el PDF (pagina 8). Son datos INFERIDOS de una
figura, con una resolucion de ~0.002 p.u.; ver datos_heydari.py para las
comprobaciones cruzadas con la Tabla 4 del paper.

Uso:   python digitalizar_fig2_heydari.py  [ruta_al_pdf]
Requiere pymupdf y numpy.
"""
import sys
from pathlib import Path

import numpy as np
import pymupdf

PDF_DEF = Path(__file__).resolve().parent.parent / "1-s2.0-S0378779624000713-main.pdf"

# Calibracion medida sobre la imagen de 900 x 1296 px (pixeles):
X0, PX_POR_H = 139.5, 0.08675          # eje x: 0 h en x=139.5; 2000 h cada 173.5 px
HB = [0, 1000, 4000, 7000, 8000, 8760]  # fronteras de los bloques (lineas punteadas)
# panel superior (viento): y=309.5 -> 0 pu ; y=1.0 -> 0.7 pu
# panel inferior (carga) : y=510.5 -> 1.0 pu ; y=1173.5 -> 0.2 pu


def _y_viento(y):
    return (309.5 - y) / (309.5 - 1.0) * 0.7


def _y_carga(y):
    return 1.0 - (y - 510.5) / (1173.5 - 510.5) * 0.8


def _centros(filas):
    out, cur = [], [filas[0]]
    for r in filas[1:]:
        if r - cur[-1] <= 1:
            cur.append(r)
        else:
            out.append(float(np.mean(cur)))
            cur = [r]
    out.append(float(np.mean(cur)))
    return out


def imagen_fig2(pdf):
    doc = pymupdf.open(str(pdf))
    pag = doc[7]                                        # pagina 8
    mejor = max(pag.get_images(full=True), key=lambda i: i[3])   # la mas alta (900 x 1296)
    pix = pymupdf.Pixmap(doc, mejor[0])
    if pix.n - pix.alpha >= 4:
        pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
    if (pix.width, pix.height) != (900, 1296):
        raise RuntimeError(f"imagen inesperada {pix.width}x{pix.height}; la calibracion es para 900x1296")
    return np.frombuffer(pix.samples, dtype=np.uint8).reshape(
        pix.height, pix.width, pix.n)[:, :, :3].mean(axis=2)


def digitalizar(pdf=PDF_DEF):
    im = imagen_fig2(pdf)
    xb = [X0 + h * PX_POR_H for h in HB]
    carga, viento = [], []
    for b in range(5):
        xa, xe = int(xb[b] + 8), int(xb[b + 1] - 8)         # dentro del bloque
        for dest, (y0, y1), cal in ((viento, (3, 308), _y_viento), (carga, (513, 1172), _y_carga)):
            frac = (im[y0:y1, xa:xe] < 80).mean(axis=1)       # segmentos negros (la curva es gris)
            filas = [y0 + i for i, f in enumerate(frac) if f > 0.85]
            dest.append([round(cal(c), 4) for c in _centros(filas)])
    return [c[0] for c in carga], viento


if __name__ == "__main__":
    pdf = Path(sys.argv[1]) if len(sys.argv) > 1 else PDF_DEF
    carga, viento = digitalizar(pdf)
    for b in range(5):
        print(f"bloque {b + 1} ({HB[b]}-{HB[b + 1]} h, {HB[b + 1] - HB[b]} h): "
              f"carga={carga[b]}  viento={viento[b]}")
    try:
        from datos_heydari import CARGA, VIENTO
        ok = (np.allclose(carga, CARGA, atol=1e-4) and np.allclose(viento, VIENTO, atol=1e-4))
        print("coincide con datos_heydari.py:", ok)
    except ImportError:
        pass
