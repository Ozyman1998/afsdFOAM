"""
Mallado con gmsh por extrusión de la planta de CadQuery (como la tesis:
malla 2D de la superficie del sustrato extruida hacia arriba y hacia abajo).

  - Bajo la herramienta: O-grid transfinito de cuadriláteros -> hexaedros
  - Fuera: triángulos -> prismas (o cuadriláteros si RECOMBINE_OUTER)
  - Capas en z graduadas, conformes entre las tres regiones

Salida: afsd.msh (formato 2.2 ASCII, en metros) con
  volúmenes físicos -> cellZones  build, substrate, anvil
  superficies físicas -> patches
"""
import math
import sys

import gmsh

from params import *

TOL = 1e-6


# ---------------------------------------------------------------- utilidades
def graded_layers(n, total, dz1):
    """Alturas acumuladas normalizadas (0..1] de n capas geométricas,
    primera capa ~dz1. Devuelve (numElements, heights) para occ.extrude."""
    if abs(n * dz1 - total) < 1e-9:
        r = 1.0
    else:
        lo, hi = 1.0 + 1e-9, 2.0   # razón r tal que dz1*(r^n-1)/(r-1) = total
        f = lambda r: dz1 * (r**n - 1) / (r - 1) - total
        if f(lo) > 0:              # dz1 demasiado grande: capas uniformes
            r = 1.0
        else:
            while f(hi) < 0:
                hi *= 2
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
            r = 0.5 * (lo + hi)
    dz = [dz1 * r**i for i in range(n)] if r != 1.0 else [total / n] * n
    s = sum(dz)
    cum, acc = [], 0.0
    for d in dz:
        acc += d
        cum.append(acc / s)
    cum[-1] = 1.0
    return [1] * n, cum


def extrude(faces, dz, n, dz1=None):
    """Extruye caras 2D en z. Devuelve (volúmenes, caras superiores)."""
    if dz1 is None:
        num, h = [n], []
    else:
        num, h = graded_layers(n, abs(dz), dz1)
    out = gmsh.model.occ.extrude(
        [(2, f) for f in faces], 0, 0, dz, numElements=num, heights=h, recombine=True
    )
    vols, tops = [], []
    for i, (d, t) in enumerate(out):
        if d == 3:
            vols.append(t)
            tops.append(out[i - 1][1])
    return vols, tops


def endpoints_r(curve):
    pts = gmsh.model.getBoundary([(1, curve)], oriented=False)
    rs = []
    for _, p in pts:
        x, y, _ = gmsh.model.getValue(0, p, [])
        rs.append(math.hypot(x, y))
    return rs


def near(a, b, tol=1e-3):
    return abs(a - b) < tol


# ---------------------------------------------------------------- planta
def build_planform():
    ents = gmsh.model.occ.importShapes(BREP_FILE)
    faces = [t for d, t in ents if d == 2]
    out, _ = gmsh.model.occ.fragment([(2, f) for f in faces], [])
    gmsh.model.occ.synchronize()
    faces = [t for d, t in out if d == 2]

    tool_faces, inlet_faces, outer_faces = [], [], []
    for f in faces:
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(2, f)
        x, y, _ = gmsh.model.occ.getCenterOfMass(2, f)
        r = math.hypot(x, y)
        if xmax > R_TOOL + 1e-3:         # rectángulo con agujero: su centroide
            outer_faces.append(f)       # cae en el origen, se distingue por bbox
        elif r < R_INLET:
            inlet_faces.append(f)
        elif r < R_TOOL:
            tool_faces.append(f)
        else:
            raise RuntimeError(f"cara de planta sin clasificar {f}")
    assert len(inlet_faces) == 5 and len(tool_faces) == 4 and len(outer_faces) == 1, (
        len(inlet_faces), len(tool_faces), len(outer_faces))
    return inlet_faces, tool_faces, outer_faces


def set_planform_mesh(inlet_faces, tool_faces, outer_faces):
    m = gmsh.model.mesh
    o_grid = inlet_faces + tool_faces
    curves = {c for f in o_grid for _, c in gmsh.model.getBoundary([(2, f)])}
    curves = {abs(c) for c in curves}

    for c in curves:
        r0, r1 = sorted(endpoints_r(c))
        if near(r0, r1):                       # arcos y lados del cuadrado
            m.setTransfiniteCurve(c, N_CIRC + 1)
        elif near(r0, R_SQUARE_CORNER) and near(r1, R_INLET):
            m.setTransfiniteCurve(c, N_RAD_INLET + 1)
        elif near(r0, R_INLET) and near(r1, R_TOOL):
            # algo más fino junto al borde exterior de la herramienta
            m.setTransfiniteCurve(c, N_RAD_TOOL + 1, "Bump", 0.5)
        else:
            raise RuntimeError(f"curva no clasificada {c}: r = {r0:.3f}, {r1:.3f}")

    for f in o_grid:
        m.setTransfiniteSurface(f)
        m.setRecombine(2, f)

    if RECOMBINE_OUTER:
        for f in outer_faces:
            m.setRecombine(2, f)

    # Tamaño fuera de la herramienta: banda fina en la trayectoria
    fld = gmsh.model.mesh.field
    box = fld.add("Box")
    fld.setNumber(box, "VIn", SIZE_PATH)
    fld.setNumber(box, "VOut", SIZE_FAR)
    fld.setNumber(box, "XMin", -L_DOMAIN)
    fld.setNumber(box, "XMax", L_DOMAIN)
    fld.setNumber(box, "YMin", -Y_PATH)
    fld.setNumber(box, "YMax", Y_PATH)
    fld.setNumber(box, "ZMin", -1e3)
    fld.setNumber(box, "ZMax", 1e3)
    fld.setNumber(box, "Thickness", 10.0)
    fld.setAsBackgroundMesh(box)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


# ---------------------------------------------------------------- volúmenes
def build_volumes(inlet_faces, tool_faces, outer_faces):
    base = inlet_faces + tool_faces + outer_faces

    gap_vols, gap_tops = extrude(base, H_LAYER, N_GAP)
    outer_top = [t for f, t in zip(base, gap_tops) if f in outer_faces]
    up_vols, _ = extrude(outer_top, H_BUILD - H_LAYER, N_UP, DZ1_UP)

    sub_vols, sub_bots = extrude(base, -T_SUB, N_SUB, DZ1_SUB)
    anv_vols, _ = extrude(sub_bots, -T_ANVIL, N_ANVIL, DZ1_ANVIL)

    # Las extrusiones de caras vecinas duplican las caras laterales:
    # fusionarlas para que la malla sea conforme.
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()


# ---------------------------------------------------------------- grupos
def classify_and_tag():
    vols = {"build": [], "substrate": [], "anvil": []}
    for _, v in gmsh.model.getEntities(3):
        _, _, z = gmsh.model.occ.getCenterOfMass(3, v)
        key = "build" if z > 0 else ("substrate" if z > -T_SUB else "anvil")
        vols[key].append(v)

    boundary = gmsh.model.getBoundary(
        [(3, v) for vs in vols.values() for v in vs], combined=True, oriented=False
    )
    region_of = lambda z: "build" if z > 0 else ("substrate" if z > -T_SUB else "anvil")

    patches = {}
    for _, s in boundary:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)
        x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
        reg = region_of(z)
        if near(zmin, zmax) and near(zmin, H_BUILD):
            name = "build_top"
        elif near(zmin, zmax) and near(zmin, -T_SUB - T_ANVIL):
            name = "anvil_bottom"
        elif near(xmin, xmax) and near(xmin, L_DOMAIN / 2):
            name = f"{reg}_front"     # x = +L/2: entra el material (marco de la herramienta)
        elif near(xmin, xmax) and near(xmin, -L_DOMAIN / 2):
            name = f"{reg}_back"      # x = -L/2: sale el depósito
        elif near(ymin, ymax) and near(abs(ymin), W_DOMAIN / 2):
            name = f"{reg}_sides"
        elif near(zmin, zmax) and near(zmin, H_LAYER):
            name = "inlet" if math.hypot(x, y) < R_INLET else "toolFace"
        elif zmin >= H_LAYER - TOL and gmsh.model.getType(2, s) == "Cylinder":
            name = "toolSide"
        else:
            raise RuntimeError(f"superficie de contorno sin clasificar {s}: "
                               f"bbox z [{zmin:.3f}, {zmax:.3f}] centro ({x:.2f},{y:.2f},{z:.2f})")
        patches.setdefault(name, []).append(s)

    for i, (name, vs) in enumerate(vols.items()):
        gmsh.model.addPhysicalGroup(3, vs, tag=1001 + i, name=name)
    for i, (name, ss) in enumerate(sorted(patches.items())):
        gmsh.model.addPhysicalGroup(2, ss, tag=1 + i, name=name)
    return vols, patches


# ---------------------------------------------------------------- main
def main():
    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("afsd")

    inlet_f, tool_f, outer_f = build_planform()
    set_planform_mesh(inlet_f, tool_f, outer_f)
    build_volumes(inlet_f, tool_f, outer_f)
    vols, patches = classify_and_tag()

    gmsh.option.setNumber("Mesh.Algorithm", 6)          # Frontal-Delaunay 2D
    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.ScalingFactor", 1e-3)   # mm -> m
    gmsh.write(MESH_FILE)

    # Resumen
    print("\n=== Resumen ===")
    for name, vs in vols.items():
        n = sum(len(gmsh.model.mesh.getElements(3, v)[1][i])
                for v in vs for i in range(len(gmsh.model.mesh.getElements(3, v)[1])))
        print(f"  {name:10s} {len(vs):3d} volúmenes  {n:8d} celdas")
    for name, ss in sorted(patches.items()):
        print(f"  patch {name:18s} {len(ss):3d} superficies")
    print(f"Escrito {MESH_FILE}")

    if "-gui" in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()


if __name__ == "__main__":
    main()
