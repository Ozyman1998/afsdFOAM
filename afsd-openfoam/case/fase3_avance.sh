#!/bin/sh
# Fase 3: avance (35-58.77 s): 50.4 mm a Vt = 2.12 mm/s. Continúa en paralelo
# desde el último tiempo de la fase 2.
set -e
cd "${0%/*}"
NP=${1:-10}

ls -d processor0/35 > /dev/null || { echo "Falta processor0/35: ¿terminó la fase 2?"; exit 1; }

foamDictionary constant/frameProperties -entry Vt -set '(0.00212 0 0)'
foamDictionary system/controlDict -entry startFrom -set latestTime
foamDictionary system/controlDict -entry endTime -set 58.77

nohup mpirun -np $NP afsdFoam -parallel > log.fase3_avance 2>&1 &
echo "Fase 3 lanzada con $NP procesos (log.fase3_avance)"
