#!/usr/bin/env python3
"""
Fase 1 (dwell): flujo de calor de la herramienta sobre el sustrato.

Como en la tesis (cap. 4.4.2): sólo sustrato y yunque; flujo de calor en el
anillo de la cara de la herramienta, proporcional al radio, con potencia total
Q = 2.05 kW. Fuera del anillo, adiabático (la convección h = 10 W/m2K sobre
la cara superior es despreciable frente a 2 kW).

    q(r) = c r,   R_i <= r <= R_o,   c = 3 Q / (2 pi (R_o^3 - R_i^3))

Lee los centros de cara del patch substrate_to_build de 0/substrate/C
(postProcess -func writeCellCentres) y añade la condición
externalWallHeatFluxTemperature (mode flux) a 0/substrate/T.
"""
import math
import os
import re
import sys

CASE = os.path.dirname(os.path.abspath(__file__))
PATCH = "substrate_to_build"

Q_TOTAL = 2050.0                                  # [W] (cap. 4)
R_O = 38.1e-3 / 2                                 # [m] radio de la herramienta
R_I = math.sqrt(9.35e-3**2 / math.pi)             # [m] radio de la entrada


def face_centres():
    txt = open(os.path.join(CASE, "0", "substrate", "C")).read()
    m = re.search(
        PATCH + r"\s*\{[^}]*?value\s+nonuniform\s+List<vector>\s*(\d+)\s*\(",
        txt,
        re.S,
    )
    if not m:
        sys.exit(f"No encuentro los centros de cara de {PATCH} en 0/substrate/C")
    n = int(m.group(1))
    vecs = re.findall(
        r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)",
        txt[m.end():],
    )[:n]
    assert len(vecs) == n, (len(vecs), n)
    return [(float(x), float(y)) for x, y, _ in vecs]


def main():
    c = 3 * Q_TOTAL / (2 * math.pi * (R_O**3 - R_I**3))
    q = []
    for x, y in face_centres():
        r = math.hypot(x, y)
        q.append(c * r if R_I <= r <= R_O else 0.0)

    block = (
        f"    // Fase 1 (dwell): flujo de la herramienta, Q = {Q_TOTAL:.0f} W\n"
        f"    {PATCH}\n    {{\n"
        "        type            externalWallHeatFluxTemperature;\n"
        "        mode            flux;\n"
        "        kappaMethod     solidThermo;\n"
        f"        q               nonuniform List<scalar> {len(q)}\n        (\n"
        + "".join(f"            {v:.8g}\n" for v in q)
        + "        );\n"
        "        value           uniform 300;\n"
        "    }\n\n"
    )

    path = os.path.join(CASE, "0", "substrate", "T")
    txt = open(path).read()
    if re.search(r"^\s*" + PATCH + r"\s*$", txt, re.M):
        sys.exit(f"0/substrate/T ya tiene una entrada {PATCH}")
    key = '    "substrate_to_.*"'
    if key not in txt:
        sys.exit("No encuentro la entrada \"substrate_to_.*\" en 0/substrate/T")
    open(path, "w").write(txt.replace(key, block + key, 1))

    n_on = sum(v > 0 for v in q)
    print(f"R_i = {R_I*1e3:.3f} mm, R_o = {R_O*1e3:.3f} mm, c = {c:.4e} W/m3")
    print(f"q_max = {c*R_O:.4e} W/m2 en r = R_o; {n_on} de {len(q)} caras con flujo")
    print(f"Escrita la condición de {PATCH} en 0/substrate/T")


if __name__ == "__main__":
    main()
