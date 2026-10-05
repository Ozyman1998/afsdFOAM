#!/usr/bin/env python3
"""
Termopares (tesis Fig. 2.5 y 2.7).

Uso:
    python3 termopares.py --sondas     escribe system/thermocouples (antes de lanzar)
    python3 termopares.py              historia de T en los termopares

Los termopares están fijos en el sustrato: centro (Mid) y a 12.7 mm a cada lado
(Adv, Rtr). En el marco de la herramienta se mueven con el sustrato:
    x_tc(t) = X0                      t <= T_ADV
    x_tc(t) = X0 - Vt (t - T_ADV)     t >  T_ADV
con X0 = 25.2 mm (mitad de los 50.4 mm depositados): por delante de la
herramienta al empezar el avance, debajo a los 12 s y detrás a los 24 s.

Lado de avance (Adv): donde la velocidad de la herramienta (ec. 13, giro
antihorario) va en el sentido del avance (+x): y < 0.

Se muestrea una red de sondas en x (cada 1 mm) y se interpola en x_tc(t).
Profundidad del termopar: no la da la tesis -> se dan tres profundidades.
"""
import glob
import os
import re
import sys

import numpy as np

CASE = os.path.dirname(os.path.abspath(__file__))

X0 = 25.2e-3          # [m]
VT = 2.12e-3          # [m/s]
T_ADV = 35.0          # [s] inicio del avance (25 s dwell + 10 s fill)
Y_TC = {"Mid": 0.0, "Adv": -12.7e-3, "Rtr": 12.7e-3}
DEPTHS = [-1.0e-3, -3.175e-3, -5.35e-3]     # [m] bajo la cara superior
XS = np.arange(-30, 31) * 1e-3              # [m]


def points():
    return [(x, y, z) for z in DEPTHS for y in Y_TC.values() for x in XS]


def write_probes():
    pts = "\n".join(f"        ({x:.6g} {y:.6g} {z:.6g})" for x, y, z in points())
    txt = f"""// Generado por termopares.py --sondas
thermocouples
{{
    type            probes;
    libs            (sampling);
    region          substrate;
    fields          (T);
    interpolationScheme cellPoint;
    writeControl    adjustableRunTime;
    writeInterval   0.1;
    probeLocations
    (
{pts}
    );
}}
"""
    open(os.path.join(CASE, "system", "thermocouples"), "w").write(txt)
    print(f"system/thermocouples: {len(points())} sondas")


def read_probes():
    files = sorted(
        glob.glob(os.path.join(CASE, "postProcessing", "**", "thermocouples", "**", "T"), recursive=True),
        key=lambda f: float(os.path.basename(os.path.dirname(f))),
    )
    if not files:
        sys.exit("No hay datos en postProcessing/substrate/thermocouples")
    data = {}
    for f in files:            # los reinicios posteriores sobrescriben
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            v = np.array(line.split(), dtype=float)
            data[v[0]] = v[1:]
    t = np.array(sorted(data))
    T = np.array([data[k] for k in t])
    return t, T


def main():
    if "--sondas" in sys.argv:
        write_probes()
        return

    t, T = read_probes()
    nx, ny = len(XS), len(Y_TC)
    T = T.reshape(len(t), len(DEPTHS), ny, nx)
    xtc = np.where(t <= T_ADV, X0, X0 - VT * (t - T_ADV))

    adv = t >= T_ADV
    for iz, z in enumerate(DEPTHS):
        out = os.path.join(CASE, f"termopares_z{abs(z)*1e3:.2f}mm.csv")
        with open(out, "w") as fo:
            fo.write("t_avance_s," + ",".join(f"{k}_C" for k in Y_TC) + "\n")
            for it in np.where(adv)[0]:
                vals = [
                    np.interp(xtc[it], XS, T[it, iz, iy]) - 273.15
                    for iy in range(ny)
                ]
                fo.write(f"{t[it]-T_ADV:.3f}," + ",".join(f"{v:.2f}" for v in vals) + "\n")
        print(f"Escrito {os.path.basename(out)}")

    # Fin del dwell y del fill: T en la posición de los termopares
    for label, tt in (("fin dwell", 25.0), ("inicio avance", T_ADV)):
        it = np.argmin(abs(t - tt))
        if abs(t[it] - tt) < 0.2:
            vals = [np.interp(X0, XS, T[it, 1, iy]) - 273.15 for iy in range(ny)]
            print(f"  {label} (t = {t[it]:.1f} s, z = -3.175 mm): "
                  + "  ".join(f"{k} {v:.0f} °C" for k, v in zip(Y_TC, vals)))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if not adv.any():
        return
    fig, axes = plt.subplots(1, len(DEPTHS), figsize=(14, 4), sharey=True)
    colors = {"Mid": "k", "Adv": "r", "Rtr": "b"}
    exp = os.path.join(CASE, "exp_optimo.csv")
    for iz, (ax, z) in enumerate(zip(axes, DEPTHS)):
        for iy, k in enumerate(Y_TC):
            Tk = [np.interp(xtc[it], XS, T[it, iz, iy]) - 273.15 for it in np.where(adv)[0]]
            ax.plot(t[adv] - T_ADV, Tk, color=colors[k], label=f"FVM {k}")
        if os.path.isfile(exp):
            d = np.genfromtxt(exp, delimiter=",", names=True, dtype=None, encoding=None)
            for k in Y_TC:
                sel = d["serie"] == k
                ax.plot(d["t"][sel], d["T"][sel], "o", ms=3, color=colors[k], label=f"Exp. {k}")
        ax.set_title(f"z = {z*1e3:.2f} mm")
        ax.set_xlabel("Tiempo de avance (s)")
        ax.set_xlim(0, 24)
        ax.set_ylim(0, 550)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("T (°C)")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(CASE, "termopares.png"), dpi=150)
    print("Escrito termopares.png")


if __name__ == "__main__":
    main()
