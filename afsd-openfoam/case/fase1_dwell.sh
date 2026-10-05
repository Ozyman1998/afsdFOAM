#!/bin/sh
# Fase 1: dwell (0-25 s). Sólo sustrato y yunque; 2.05 kW en el anillo de la
# herramienta (tesis cap. 4). En serie: son ~90k celdas sin fluido.
set -e
cd "${0%/*}"

foamListTimes -rm
rm -rf 0 postProcessing processor*
cp -r 0.orig 0

cp constant/regionProperties.dwell constant/regionProperties
foamDictionary constant/frameProperties -entry Vt -set '(0 0 0)'

python3 termopares.py --sondas

postProcess -region substrate -func writeCellCentres -time 0 > log.writeCellCentres 2>&1
python3 dwell_flux.py
rm -f 0/substrate/C 0/substrate/Cx 0/substrate/Cy 0/substrate/Cz

foamDictionary system/controlDict -entry startFrom -set startTime
foamDictionary system/controlDict -entry endTime -set 25
foamDictionary system/controlDict -entry maxDeltaT -set 0.05
foamDictionary system/controlDict -entry functions/metalVolume/enabled -set false

afsdFoam > log.fase1_dwell 2>&1
tail -3 log.fase1_dwell
