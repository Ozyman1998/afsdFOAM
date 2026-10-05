"""
Comprueba la malla antes de pasarla a OpenFOAM:
  - tipos de celda (solo hexaedros y prismas)
  - nodos duplicados (malla no conforme entre volúmenes/regiones)
  - toda cara de contorno pertenece a un patch, y ningún patch es interno
"""
import sys
from collections import Counter

import gmsh
import numpy as np

from params import MESH_FILE

# caras de cada tipo de celda gmsh (índices locales de nodos)
FACES = {
    5: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    6: [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    4: [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)],
    7: [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
NAMES = {4: "tetra", 5: "hexa", 6: "prisma", 7: "pirámide"}

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(MESH_FILE)
ok = True

# nodos duplicados
tags, xyz, _ = gmsh.model.mesh.getNodes()
xyz = np.round(xyz.reshape(-1, 3) / 1e-9).astype(np.int64)
dup = len(xyz) - len(np.unique(xyz, axis=0))
print(f"Nodos: {len(tags)}  duplicados: {dup}")
ok &= dup == 0

# celdas y caras
faces = Counter()
types = Counter()
for etype, _, enodes in zip(*gmsh.model.mesh.getElements(3)):
    nn = gmsh.model.mesh.getElementProperties(etype)[3]
    conn = enodes.reshape(-1, nn)
    types[NAMES.get(etype, etype)] += len(conn)
    for loc in FACES[etype]:
        for f in conn[:, loc]:
            faces[tuple(sorted(f))] += 1
print("Celdas:", dict(types), " total", sum(types.values()))
ok &= set(types) <= {"hexa", "prisma"}

bad = sum(1 for c in faces.values() if c > 2)
bnd = {f for f, c in faces.items() if c == 1}
print(f"Caras: {len(faces)}  de contorno: {len(bnd)}  compartidas por >2 celdas: {bad}")
ok &= bad == 0

# caras etiquetadas en patches
tagged = set()
for dim, ptag in gmsh.model.getPhysicalGroups(2):
    name = gmsh.model.getPhysicalName(dim, ptag)
    n_in = n_out = 0
    for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, ptag):
        for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2, ent)):
            nn = gmsh.model.mesh.getElementProperties(etype)[3]
            for f in enodes.reshape(-1, nn):
                key = tuple(sorted(f))
                tagged.add(key)
                n_in += key in bnd
                n_out += key not in bnd
    print(f"  {name:18s} caras {n_in:7d}" + (f"   ¡{n_out} internas!" if n_out else ""))
    ok &= n_out == 0
untagged = len(bnd - tagged)
print(f"Caras de contorno sin patch: {untagged}")
ok &= untagged == 0

gmsh.finalize()
print("OK" if ok else "FALLO")
sys.exit(0 if ok else 1)
