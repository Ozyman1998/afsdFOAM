#!/bin/sh
# Geometría (CadQuery) + malla (gmsh) + paso a OpenFOAM con 3 regiones.
# Necesita el entorno de Python con cadquery y gmsh activado.
set -e
cd "${0%/*}"

echo "== Geometría y malla"
( cd geometry && python3 geometry.py && python3 mesh.py > log.mesh && python3 check_mesh.py )

echo "== OpenFOAM"
rm -rf constant/polyMesh constant/build/polyMesh constant/substrate/polyMesh constant/anvil/polyMesh
gmshToFoam geometry/afsd.msh > log.gmshToFoam 2>&1
splitMeshRegions -cellZones -overwrite > log.splitMeshRegions 2>&1
for r in build substrate anvil; do
    checkMesh -region $r > log.checkMesh.$r 2>&1
done
grep -H -E "Mesh OK|Failed" log.checkMesh.*
ls constant
