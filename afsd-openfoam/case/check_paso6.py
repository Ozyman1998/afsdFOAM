#!/usr/bin/env python3
"""
Paso 6: verifica la viscosidad de Sellars-Tegart que calcula afsdFoam.

Lee T, epsDot y mu.metal de la región build en un tiempo escrito, recalcula
mu con la ley de la tesis (ec. 6-9) y los límites de phaseProperties, y
compara celda a celda.

Uso:  python3 check_paso6.py [tiempo]      (por defecto, el último)
"""
import math
import os
import re
import sys

import numpy as np

CASE = os.path.dirname(os.path.abspath(__file__))


def strip_comments(txt):
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    return re.sub(r"//.*", "", txt)


def coeffs():
    txt = strip_comments(
        open(os.path.join(CASE, "constant", "build", "phaseProperties")).read()
    )
    m = re.search(r"SellarsTegartCoeffs\s*\{(.*?)\}", txt, re.S)
    if not m:
        sys.exit("No hay SellarsTegartCoeffs en phaseProperties")
    c = {}
    for k in ("sigmaN", "A", "n", "Q", "R", "epsDotMin", "muMax", "Tmin", "Tmax"):
        mm = re.search(rf"\b{k}\s+([-+0-9.eE]+)\s*;", m.group(1))
        if not mm:
            sys.exit(f"Falta {k} en SellarsTegartCoeffs")
        c[k] = float(mm.group(1))
    return c


def read_internal(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", txt)
    if m:
        return float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt)
    if not m:
        sys.exit(f"No se puede leer internalField de {path}")
    n = int(m.group(1))
    start = m.end()
    end = txt.index(")", start)
    vals = np.array(txt[start:end].split(), dtype=float)
    assert vals.size == n, (vals.size, n)
    return vals


def latest_time():
    times = []
    for d in os.listdir(CASE):
        try:
            t = float(d)
        except ValueError:
            continue
        if t > 0 and os.path.isfile(os.path.join(CASE, d, "build", "mu.metal")):
            times.append((t, d))
    if not times:
        sys.exit("No hay tiempos escritos con build/mu.metal")
    return max(times)[1]


def main():
    tdir = sys.argv[1] if len(sys.argv) > 1 else latest_time()
    c = coeffs()
    base = os.path.join(CASE, tdir, "build")
    T = read_internal(os.path.join(base, "T"))
    e = read_internal(os.path.join(base, "epsDot"))
    mu = read_internal(os.path.join(base, "mu.metal"))
    n = max(np.size(T), np.size(e), np.size(mu))
    T, e, mu = (np.broadcast_to(a, n).astype(float) for a in (T, e, mu))

    ee = np.maximum(e, c["epsDotMin"])
    Te = np.clip(T, c["Tmin"], c["Tmax"])
    Z = ee * np.exp(c["Q"] / (c["R"] * Te))
    sigma = c["sigmaN"] * np.arcsinh((Z / c["A"]) ** (1.0 / c["n"]))
    mu_ref = np.minimum(sigma / (3.0 * ee), c["muMax"])

    rel = np.abs(mu - mu_ref) / mu_ref
    capped = np.mean(mu_ref >= c["muMax"] * (1 - 1e-12))

    print(f"Tiempo {tdir}: {n} celdas")
    print(f"  T      : {T.min():.1f} - {T.max():.1f} K")
    print(f"  epsDot : {e.min():.3e} - {e.max():.3e} 1/s")
    print(f"  mu.metal: {mu.min():.3e} - {mu.max():.3e} Pa s"
          f"  ({100*capped:.1f} % de celdas en muMax)")
    print(f"  Error relativo frente a la ley: max {rel.max():.2e}, medio {rel.mean():.2e}")

    # Puntos de referencia de la ley (independientes de la malla)
    print("\n  Referencia: mu [Pa s] para epsDot = 1, 10, 100 1/s")
    for Tk in (600, 700, 800, 850):
        row = []
        for eps in (1.0, 10.0, 100.0):
            Zr = eps * math.exp(c["Q"] / (c["R"] * Tk))
            s = c["sigmaN"] * math.asinh((Zr / c["A"]) ** (1 / c["n"]))
            row.append(f"{min(s/(3*eps), c['muMax']):9.2e}")
        print(f"    T = {Tk} K: " + "  ".join(row))

    ok = rel.max() < 1e-5   # campos escritos con writePrecision 8
    print("\nPASS" if ok else "\nFAIL (error > 1e-5)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
