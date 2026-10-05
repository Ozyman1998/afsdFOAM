#!/bin/sh
# Fase 2: fill (25-35 s). Modelo completo; la herramienta gira y alimenta
# sin avanzar. En paralelo.
set -e
cd "${0%/*}"
NP=${1:-10}
T0=25

[ -d "$T0/substrate" ] || { echo "Falta $T0/substrate: ¿terminó la fase 1?"; exit 1; }

cp constant/regionProperties.full constant/regionProperties
cp -r 0.orig/build "$T0/build"

# Cara superior del sustrato: de flujo de la herramienta a acoplada con build
foamDictionary "$T0/substrate/T" -entry boundaryField/substrate_to_build -set \
  '{ type compressible::turbulentTemperatureCoupledBaffleMixed; Tnbr T; kappaMethod solidThermo; value uniform 300; }'

# El dwell terminó con pasos grandes: arrancar el fluido con uno pequeño
foamDictionary "$T0/uniform/time" -entry deltaT -set 1e-4
foamDictionary "$T0/uniform/time" -entry deltaT0 -set 1e-4

foamDictionary constant/frameProperties -entry Vt -set '(0 0 0)'
foamDictionary constant/build/toolProperties -entry feedVelocity -set 0.00212
foamDictionary system/controlDict -entry startFrom -set latestTime
foamDictionary system/controlDict -entry endTime -set 35
foamDictionary system/controlDict -entry maxDeltaT -set 0.005
foamDictionary system/controlDict -entry functions/metalVolume/enabled -set true

for d in system/decomposeParDict system/build/decomposeParDict \
         system/substrate/decomposeParDict system/anvil/decomposeParDict; do
  foamDictionary $d -entry numberOfSubdomains -set $NP
done
rm -rf processor*
decomposePar -allRegions -latestTime > log.decomposePar.fill 2>&1
tail -2 log.decomposePar.fill

nohup mpirun -np $NP afsdFoam -parallel > log.fase2_fill 2>&1 &
echo "Fase 2 lanzada con $NP procesos (log.fase2_fill)"
