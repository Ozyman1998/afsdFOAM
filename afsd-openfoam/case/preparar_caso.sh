#!/bin/sh
# Convierte una copia NUEVA del tutorial multiRegionHeater en el caso AFSD.
# Ejecutar desde la carpeta del caso, DESPUÉS de copiar encima case/ del paquete.
#   ./preparar_caso.sh [N]      N = núcleos físicos para el paralelo (10 por defecto)
set -e
cd "${0%/*}"
NP=${1:-10}
TUT=$FOAM_TUTORIALS/heatTransfer/chtMultiRegionFoam/multiRegionHeater

echo "== 1. Sólidos: propiedades y esquemas del 'heater' del tutorial"
for r in substrate anvil; do
    mkdir -p constant/$r system/$r
    find $TUT/constant/heater -maxdepth 1 \( -type f -o -type l \) -exec cp -L {} constant/$r/ \;
    cp -L $TUT/system/heater/fvSchemes $TUT/system/heater/fvSolution system/$r/
done

echo "== 2. Quitar lo que es sólo del tutorial"
rm -rf constant/topAir constant/bottomWater constant/heater constant/leftSolid \
       constant/rightSolid constant/cellToRegion constant/polyMesh
rm -rf system/topAir system/bottomWater system/heater system/leftSolid system/rightSolid
rm -f  system/vtkWrite system/README system/blockMeshDict system/topoSetDict
rm -f  0.orig/epsilon 0.orig/k 0.orig/p 0.orig/p_rgh 0.orig/T 0.orig/U
rm -f  Allrun Allrun.pre Allclean

echo "== 3. Propiedades de los sólidos (constantes)"
#   Sustrato AA6061-T6
foamDictionary constant/substrate/thermophysicalProperties -entry mixture/transport/kappa     -set 167
foamDictionary constant/substrate/thermophysicalProperties -entry mixture/thermodynamics/Cp   -set 896
foamDictionary constant/substrate/thermophysicalProperties -entry mixture/equationOfState/rho -set 2700
#   Yunque, acero al carbono
foamDictionary constant/anvil/thermophysicalProperties     -entry mixture/transport/kappa     -set 55
foamDictionary constant/anvil/thermophysicalProperties     -entry mixture/thermodynamics/Cp   -set 486
foamDictionary constant/anvil/thermophysicalProperties     -entry mixture/equationOfState/rho -set 7850

echo "== 4. Sólidos: esquema de advección y solver no simétrico"
for r in substrate anvil; do
    foamDictionary system/$r/fvSchemes -entry 'divSchemes/div(phiSolid,h)' -set 'Gauss linearUpwind grad(h)'
    for f in h hFinal; do
        foamDictionary system/$r/fvSolution -entry solvers/$f/solver         -set PBiCGStab
        foamDictionary system/$r/fvSolution -entry solvers/$f/preconditioner -set DILU
    done
done

echo "== 5. Bucles extra de energía (tesis 2.2.3; cap. 5: sin beneficio -> 0)"
foamDictionary system/fvSolution -entry PIMPLE/nEnergyCouplingLoops -set 0

echo "== 6. Paralelo: $NP subdominios en todas las regiones"
foamDictionary system/decomposeParDict -entry numberOfSubdomains -set $NP
foamDictionary system/decomposeParDict -entry method -set scotch
for r in build substrate anvil; do cp system/decomposeParDict system/$r/; done

echo "== 7. Sondas de los termopares"
python3 termopares.py --sondas

chmod +x *.sh *.py geometry/*.py
echo "Caso preparado. Siguiente: ./mallar.sh"
