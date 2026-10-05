# afsd-openfoam

Réplica en **OpenFOAM v2412 (ESI)** del modelo de volúmenes finitos de *additive friction stir deposition* (AFSD) del capítulo 2 de:

> K. C. Kincaid, *Process Modeling for Additive Friction Stir Deposition via the Finite Volume Method for Terrestrial and Non-Terrestrial Applications*, PhD dissertation, The University of Alabama, 2022.

Incluye un solver propio (`afsdFoam`), una geometría y una malla paramétricas (CadQuery + gmsh) y la secuencia completa del proceso, **dwell → fill → avance**, con la comparación con los termopares de la tesis.

> **Estado:** el modelo está validado por partes (ver [Estado](#estado)). La simulación completa del caso óptimo y la comparación con la Fig. 2.7 de la tesis están en curso.

---

## Contenido

- [El modelo](#el-modelo)
- [Estado](#estado)
- [Antes de empezar](#antes-de-empezar)
- [Paso 1 — Compilar el solver](#paso-1--compilar-el-solver)
- [Paso 2 — Preparar el caso](#paso-2--preparar-el-caso)
- [Paso 3 — Geometría, malla y regiones](#paso-3--geometría-malla-y-regiones)
- [Paso 4 — Fase 1: dwell](#paso-4--fase-1-dwell-0--25-s)
- [Paso 5 — Fase 2: fill](#paso-5--fase-2-fill-25--35-s)
- [Paso 6 — Fase 3: avance](#paso-6--fase-3-avance-35--5877-s)
- [Paso 7 — Resultados](#paso-7--resultados)
- [Parar y reanudar un cálculo](#parar-y-reanudar-un-cálculo)
- [Qué hace cada fichero](#qué-hace-cada-fichero)
- [Decisiones y supuestos frente a la tesis](#decisiones-y-supuestos-frente-a-la-tesis)
- [Problemas frecuentes](#problemas-frecuentes)
- [Pendiente](#pendiente)
- [Licencia](#licencia)

---

## El modelo

Tres regiones resueltas en el marco de la herramienta, que avanza a Vt = 2.12 mm/s. El material de aporte y el aire forman un fluido VOF; el sustrato y el yunque son sólidos que se trasladan a −Vt.

| Región | Qué se resuelve | Ecuaciones de la tesis |
| --- | --- | --- |
| `build` (metal AA6061 + aire, 0–50 mm) | VOF incompresible, cantidad de movimiento con viscosidad de Sellars–Tegart, energía con calentamiento viscoso | ec. 1–11 |
| `substrate` (AA6061, 152.4 × 76.2 × 6.35 mm) | Conducción + advección a −Vt | ec. 19 |
| `anvil` (acero, 12.7 mm) | Conducción + advección a −Vt | ec. 19 |
| Herramienta (Ø 38.1 mm, a 1.5 mm del sustrato) | Condiciones de contorno: giro a 300 rpm con deslizamiento δ = 0.25, entrada de la barra, calor de fricción, conducción por la barra | ec. 12–17 |

Algoritmo por paso de tiempo (tesis 2.2.3): PIMPLE en `build` hasta tolerancia → energía en los sólidos → bucles extra solo de energía (desactivados por defecto).

## Estado

| Parte | Estado | Cómo se comprobó |
| --- | --- | --- |
| Solver base, geometría, malla y 3 regiones | Validado | `check_mesh.py` OK; `checkMesh` OK en las tres regiones |
| Advección en sólidos | Validado | Un bloque caliente se desplaza exactamente Vt·t |
| VOF metal/aire | Validado | Volumen de metal = caudal × t (error 0.00 %), α en [0, 1] |
| Viscosidad de Sellars–Tegart | Validado | `check_paso6.py`: error 7·10⁻⁸ frente a la ley recalculada |
| Herramienta | Validado | Estable 0.5 s; giro y llenado comprobados en ParaView |
| Fase 1 (dwell) | Ejecutada | Termopares al final: Mid 152 °C, Adv 143 °C, Rtr 138 °C |
| Fase 2 (fill) | En curso | — |
| Fase 3 (avance) y comparación con la Fig. 2.7 | Pendiente | — |

---

## Antes de empezar

Necesitas OpenFOAM **v2412 de ESI** (openfoam.com, no la versión .org) y un entorno de Python para la geometría.

**Requisitos:**

- OpenFOAM v2412 cargado en la terminal (el prompt empieza por `openfoam2412:`).
- Python con CadQuery, gmsh, numpy y matplotlib. Instalación, una sola vez:

```bash
sudo apt install libglu1-mesa libxcursor1 libxinerama1 libxft2
python3 -m venv ~/venvs/afsd
source ~/venvs/afsd/bin/activate
pip install cadquery gmsh numpy matplotlib
```

- Opcional: Spyder en el mismo entorno (`pip install spyder`) para ver la geometría, y ParaView para los resultados.

**Descargar el repositorio:**

```bash
git clone https://github.com/<usuario>/afsd-openfoam.git ~/afsd-openfoam
```

**Reglas de oro.** Todas salen de errores reales:

1. **Copia los comandos desde el ordenador, no desde el móvil.** Las comillas tipográficas (“ ”) rompen los ficheros de OpenFOAM.
2. **Un comando por línea.** Si pegas dos juntos, el segundo se convierte en argumentos del primero.
3. **Para borrar resultados, usa `foamListTimes -rm`, nunca `rm -rf [0-9]*`.** Ese patrón también borra `0.orig`.
4. **Lanza los cálculos largos con `nohup … &` y comprueba con `pgrep -a afsdFoam` que hay uno solo.** Dos procesos en el mismo caso estropean los resultados y el log.
5. **Usa `grep -a`** para leer los logs: a veces contienen bytes raros y `grep` los trata como binarios.
6. **Con números decimales en `awk`, antepón `LC_ALL=C`.** Con la configuración regional española, `awk` usa la coma decimal.

---

## Paso 1 — Compilar el solver

`afsdFoam` parte de una copia de `chtMultiRegionFoam`. Se le borra la parte fluida (gas compresible) y se le pone encima la de este repositorio. Las cabeceras de los sólidos y del bucle multirregión se reutilizan del original.

```bash
mkdir -p $WM_PROJECT_USER_DIR/applications/solvers
cd $WM_PROJECT_USER_DIR/applications/solvers
cp -r $FOAM_SOLVERS/heatTransfer/chtMultiRegionFoam afsdFoam
cd afsdFoam
rm -rf chtMultiRegionSimpleFoam chtMultiRegionTwoPhaseEulerFoam fluid
rm -f chtMultiRegionFoam.C
cp -r ~/afsd-openfoam/solver/. .
wclean
wmake > log.wmake 2>&1; tail -2 log.wmake; grep -n -m 10 " error" log.wmake
ls $FOAM_USER_APPBIN
```

**Cómo saber que ha ido bien:** la última línea de `tail` termina en `.../bin/afsdFoam`, el `grep` no muestra nada y `ls` muestra `afsdFoam`. Los *warnings* de variables sin usar son normales.

## Paso 2 — Preparar el caso

El caso parte del tutorial `multiRegionHeater`, porque de él salen las propiedades termofísicas y los esquemas numéricos de los sólidos. Encima se copian los ficheros de `case/`, y `preparar_caso.sh` hace el resto. Usa una carpeta nueva.

```bash
mkdir -p $WM_PROJECT_USER_DIR/run && cd $WM_PROJECT_USER_DIR/run
cp -r $FOAM_TUTORIALS/heatTransfer/chtMultiRegionFoam/multiRegionHeater afsdCase
cd afsdCase
cp -r ~/afsd-openfoam/case/. .
./preparar_caso.sh 10
```

El `10` es el número de procesos para el paralelo. Debe ser el número de **núcleos físicos**, no de hilos:

```bash
LC_ALL=C lscpu | grep -E "^Core\(s\) per socket|^Socket\(s\)"
```

Núcleos físicos = *Core(s) per socket* × *Socket(s)*.

**Qué hace `preparar_caso.sh`:**

1. Copia al sustrato y al yunque las propiedades y los esquemas del `heater` del tutorial.
2. Borra las regiones y ficheros del tutorial que no se usan.
3. Fija las propiedades de los sólidos: sustrato k = 167 W/m·K, cp = 896 J/kg·K, ρ = 2700 kg/m³; yunque k = 55, cp = 486, ρ = 7850.
4. Añade en los sólidos el esquema de la advección y un solver lineal para matrices no simétricas (PBiCGStab).
5. Desactiva los bucles extra de energía (`nEnergyCouplingLoops 0`).
6. Configura el paralelo en las tres regiones.
7. Genera las 549 sondas de los termopares (`system/thermocouples`).

**Cómo saber que ha ido bien:** termina en `Caso preparado`. Si alguna orden falla, el script se para ahí y el mensaje dice cuál.

## Paso 3 — Geometría, malla y regiones

La geometría sale de CadQuery y la malla de gmsh, que extruye la planta en vertical como en la tesis: hexaedros bajo la herramienta y prismas fuera. Las dimensiones y la densidad de malla están en `geometry/params.py`. Resultado: 232 410 celdas (build 142 266, sustrato 56 340, yunque 33 804).

```bash
source ~/venvs/afsd/bin/activate
./mallar.sh
```

**Cómo saber que ha ido bien:**

- `check_mesh.py` termina en `OK`.
- El `grep` final muestra `Mesh OK` en las tres regiones.
- `ls constant` muestra las carpetas `build`, `substrate` y `anvil`.

**Ver la geometría antes de mallar (opcional, en Spyder):** abre `geometry/geometry.py`, pulsa F5 y en la consola ejecuta `from cadquery.vis import show; show(domain_3d(), edges=True)`. Para ver la malla, ejecuta `mesh.py` con el argumento `-gui`.

## Paso 4 — Fase 1: dwell (0 → 25 s)

La herramienta gira apoyada sobre el sustrato, sin alimentar ni avanzar, y lo precalienta. Como en el cap. 4 de la tesis, solo se resuelven el sustrato y el yunque, con 2.05 kW repartidos en el anillo de la herramienta (q ∝ r). Corre en serie y tarda unos minutos.

```bash
./fase1_dwell.sh
python3 termopares.py
```

El script no imprime nada mientras calcula. Para ver por dónde va, en otra terminal: `grep -a "^Time = " log.fase1_dwell | tail -1`.

**Cómo saber que ha ido bien:** `log.fase1_dwell` termina en `End`, existe la carpeta `25` y `termopares.py` imprime algo como `fin dwell (t = 25.0 s, z = -3.175 mm): Mid 152 °C  Adv 143 °C  Rtr 138 °C`. Los dos avisos sobre `turbulentTemperatureCoupledBaffleMixed` son inofensivos.

## Paso 5 — Fase 2: fill (25 → 35 s)

La herramienta ya está a 1.5 mm del sustrato: gira y alimenta metal a 2.12 mm/s sin avanzar, hasta llenar el hueco bajo la cara. Modelo completo, en paralelo. Con 10 núcleos tarda unas **5 horas**.

```bash
./fase2_fill.sh 10
```

Puedes cerrar la terminal; no apagues ni suspendas el ordenador. **Vigilancia:**

```bash
grep -a "^Time = " log.fase2_fill | tail -1
grep -a "converged in" log.fase2_fill | tail -1
grep -a "Tool heat" log.fase2_fill | tail -1
awk '/Solving for fluid region/{f=1} /Solving for solid region/{f=0} f && /Min\/max T/{l=$0} END{print l}' log.fase2_fill
```

| Línea | Qué debe verse |
| --- | --- |
| `Time` | Sube de 25 a 35 |
| `converged in` | 2 iteraciones externas |
| `Tool heat` | Crece de ~0 W hacia el orden de kW a medida que el metal llena el hueco |
| `Min/max T` de `build` | Por debajo del solidus, 855 K |

**Al terminar** (`pgrep -a afsdFoam` vacío y el log acabado en `End`), `Tool heat` debe estar en kW. Si se queda en decenas de W, el hueco no se ha llenado: alarga el fill antes de avanzar.

## Paso 6 — Fase 3: avance (35 → 58.77 s)

La herramienta avanza 50.4 mm a 2.12 mm/s alimentando metal (Tabla 2.2, caso óptimo). Continúa desde el último tiempo de la fase 2. Tarda unas **12 horas**.

```bash
./fase3_avance.sh 10
```

Vigílala igual que la fase 2, con `log.fase3_avance`. `Time` debe llegar a 58.77.

## Paso 7 — Resultados

La validación de la tesis es la Fig. 2.7(e): la temperatura en tres termopares del sustrato durante el avance. Mid está en el centro; Adv, a 12.7 mm en el lado de avance (y < 0); Rtr, a 12.7 mm en el lado de retroceso (y > 0).

```bash
python3 termopares.py
```

Escribe `termopares_z1.00mm.csv`, `termopares_z3.17mm.csv` y `termopares_z5.35mm.csv` (tiempo de avance y T en °C de Mid, Adv y Rtr), y la gráfica `termopares.png`. Hay tres profundidades porque la tesis no da la de los termopares.

**Datos experimentales:** digitaliza la Fig. 2.7(e) con [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer) y guárdalos en `exp_optimo.csv`, en la carpeta del caso. `termopares.py` los dibuja automáticamente:

```csv
serie,t,T
Mid,0.0,350
Adv,0.0,330
Rtr,0.0,335
```

**ParaView:**

```bash
reconstructPar -allRegions > log.reconstructPar 2>&1
touch afsdCase.foam; paraview afsdCase.foam &
```

1. En *Mesh Regions* marca `build/internalMesh`, `substrate/internalMesh` y `build/patch/toolFace`, y en *Cell Arrays* los campos. Pulsa *Apply*.
2. Ve al último tiempo con *Last Frame*.
3. *Slice* vertical (Origin (0, 0, 0), Normal (0, 1, 0)) coloreado por `T` o por `alpha.metal`.
4. *Slice* horizontal a media altura del hueco (Origin (0, 0, 0.00075), Normal (0, 0, 1)) coloreado por `U`: el giro debe ser antihorario.
5. *Contour* de `alpha.metal` = 0.5 (la tesis usa 0.9): la forma del depósito.

---

## Parar y reanudar un cálculo

Los tiempos se escriben cada 1 s simulado. Para parar sin perder nada:

```bash
foamDictionary system/controlDict -entry stopAt -set writeNow
# esperar a que `pgrep -a afsdFoam` no muestre nada
foamDictionary system/controlDict -entry stopAt -set endTime
```

**Reanudar:** relanza solo el solver, que continúa desde el último tiempo escrito. **No vuelvas a ejecutar `fase2_fill.sh`**: copiaría `build` vacío otra vez en t = 25.

```bash
nohup mpirun -np 10 afsdFoam -parallel > log.fase2_fill_b 2>&1 &
```

Usa un nombre de log nuevo (`_b`, `_c`…). `termopares.py` une solo los datos de las sondas de todos los tramos.

---

## Qué hace cada fichero

**Solver** (`solver/`, sobre la copia de `chtMultiRegionFoam`):

| Fichero | Qué hace |
| --- | --- |
| `afsd.C` | Bucle de tiempo: PIMPLE del fluido hasta tolerancia → sólidos → bucles extra de energía (tesis 2.2.3) |
| `fluid/createFluidFields.H` | Campos de `build`; lee `phaseProperties` y `toolProperties` |
| `fluid/alphaEqn.H` | Fracción de metal con VOF y MULES (ec. 1) |
| `fluid/UEqn.H`, `fluid/pEqn.H` | Cantidad de movimiento y presión (ec. 3–4); interruptor `transposeStress` |
| `fluid/TEqn.H` | Energía con calentamiento viscoso 3βμε̇² (ec. 10) |
| `fluid/updateFluidProperties.H` | Mezcla de propiedades (ec. 2), Sellars–Tegart (ec. 6–9), k(T) (ec. 11) |
| `fluid/setToolVelocity.H` | U en la herramienta, la barra y la superficie del sustrato (ec. 12–13) |
| `fluid/setToolHeat.H` | Calor de fricción (ec. 15–17) y conducción por la barra (ec. 14); imprime `Tool heat` |
| `solid/solveSolid.H`, `solid/createSolidFields.H` | Conducción + advección a −Vt en los sólidos (ec. 19); lee `frameProperties` |

**Caso** (`case/`):

| Fichero | Qué hace |
| --- | --- |
| `geometry/params.py` | Todas las dimensiones y densidades de malla, con su fuente en la tesis |
| `geometry/geometry.py`, `mesh.py`, `check_mesh.py` | Geometría (CadQuery), malla (gmsh) y verificación de conformidad |
| `constant/frameProperties` | Vt, la velocidad de avance (0 en dwell y fill) |
| `constant/build/phaseProperties` | ρ, cp, k(T) de cada fase; Sellars–Tegart y sus límites; β_TQ |
| `constant/build/toolProperties` | ω, δ, Vᵢ, fricción, fracción de calor, L de la barra |
| `constant/regionProperties.full` / `.dwell` | Regiones del modelo completo / solo sólidos |
| `system/build/fvSolution` | Subrelajación, criterio del bucle externo, `transposeStress` |
| `preparar_caso.sh`, `mallar.sh` | Pasos 2 y 3 |
| `fase1_dwell.sh`, `fase2_fill.sh`, `fase3_avance.sh` | Las tres fases del proceso |
| `dwell_flux.py` | Flujo de 2.05 kW en el anillo para el dwell |
| `termopares.py` | Sondas (`--sondas`) e historia de los termopares |
| `check_paso6.py` | Comprueba que `mu.metal` sigue la ley de Sellars–Tegart |

## Decisiones y supuestos frente a la tesis

Donde la tesis no da un dato, se ha tomado uno justificado. Donde el modelo se aparta de ella, es por estabilidad numérica.

| Tema | Valor usado | Fuente o motivo | Dónde se cambia |
| --- | --- | --- | --- |
| ρ y cp del AA6061 | 2700 kg/m³, 896 J/kg·K | Supuesto (la tesis no los da) | `phaseProperties`, sólidos |
| ρ y cp del acero del yunque | 7850 kg/m³, 486 J/kg·K | Supuesto | `constant/anvil/thermophysicalProperties` |
| Yunque | 12.7 mm, acero, k = 55 | Cap. 4 (espesor), cap. 2 (k) | `params.py` |
| Altura de `build` | 50 mm | Cap. 4 | `params.py` |
| Barra cuadrada → círculo de igual área | Lado 9.35 mm → Ø 10.55 mm | Cap. 3 | `params.py` |
| Viscosidad del aire | 1.2 Pa·s (ν = 1 m²/s) | Cap. 2.3.1, artificial por estabilidad | `phaseProperties` |
| Límites de Sellars–Tegart | ε̇ ≥ 10⁻³ s⁻¹, μ ≤ 10⁷ Pa·s, T en 293–855 K | Supuesto numérico; 855 K = solidus | `phaseProperties` |
| k(T) | Evaluada con T limitada a 0–600 °C | Fuera de su rango, la cúbica da k < 0 | `phaseProperties` |
| ρcp de la mezcla | Media de (ρcp) de cada fase | Conserva la energía con el flujo VOF; la tesis promedia ρ y cp por separado | Código |
| Término μ(∇U)ᵀ | Solo dentro del metal: α·(∇U)ᵀ·∇μ_metal | Con la μ de la mezcla explota en la interfaz metal/aire | `transposeStress` |
| Presión en el calor de fricción | p local en la cara, p ≥ 0 | La tesis no dice cómo la toma | Código |
| τ en el calor de fricción | σ_f(T, 0.1 s⁻¹), sin /√3 | Cap. 2 | `toolProperties` |
| Coeficiente de fricción | μ_f = 0.4·exp(−λδωr), λ = 1 s/m | Cap. 3 (Nandan) | `toolProperties` |
| Fracción de calor al metal | 0.696 | Cap. 5 | `toolProperties` |
| L de la barra (ec. 14) | 32.5 mm = α/Vᵢ | Supuesto: longitud de penetración térmica | `toolProperties` |
| Velocidad en la entrada | Con el término −Vt de la ec. 12 | Literal de la tesis (efecto < 2 %) | Código |
| Dwell | Solo sólidos, 2.05 kW, q ∝ r | Cap. 4 | `dwell_flux.py` |
| Fill | 10 s | Cap. 4; el hueco tarda ~9 s en llenarse | `fase2_fill.sh` |
| Subrelajación | p 0.3, U 0.7, también en la iteración final | Sin ella, PIMPLE diverge con ν = 1 m²/s | `system/build/fvSolution` |
| Courant máximo | 0.25 | Cap. 3 | `controlDict` |

## Problemas frecuentes

| Síntoma | Causa | Solución |
| --- | --- | --- |
| `Invalid first character found` al leer un diccionario | Comillas tipográficas al copiar desde el móvil | Reescribe el fichero o usa `foamDictionary` |
| El solver arranca en un tiempo que no esperabas | `startFrom latestTime` y quedan tiempos viejos | `foamListTimes -rm` (nunca `rm -rf [0-9]*`) |
| `grep: coincidencia en fichero binario` | Dos `afsdFoam` escribiendo en el mismo log | `pgrep -a afsdFoam`, `pkill afsdFoam` y relanza uno solo |
| `There are not enough slots available` | Más procesos que núcleos físicos | Descompón con tantos subdominios como núcleos físicos |
| `cannot find file … system/build/decomposeParDict` | En v2412, `-allRegions` busca un diccionario por región | `preparar_caso.sh` ya los copia; si cambias N, cámbialo en los cuatro |
| `Unknown asymmetric matrix solver type PCG` | La advección hace la matriz de los sólidos no simétrica | PBiCGStab + DILU (ya en `preparar_caso.sh`) |
| `deltaT` se desploma a ~10⁻⁹ s | Inestabilidad: revisa la subrelajación y `transposeStress` | Usa los valores de `system/build/fvSolution` del repositorio |
| Cada paso agota las iteraciones externas | Tolerancia de p demasiado exigente | `outerCorrectorResidual 1e-3`, `nOuterCorrectors 10` |
| `Min/max T` parece absurdo | Los sólidos imprimen la misma etiqueta | Usa el `awk` que filtra solo la región fluida |
| El log se corta sin error | El proceso murió al cerrar la terminal | Lanza siempre con `nohup … &` |
| `awk` da coordenadas absurdas o con coma | Configuración regional española | Antepón `LC_ALL=C` |
| Avisos sobre `turbulentTemperatureCoupledBaffleMixed` | Condición con una sucesora más completa | Ignóralos; el resultado es idéntico |

## Pendiente

- [ ] Completar la fase 2 (fill) y comprobar que el calor de la herramienta llega a kW.
- [ ] Ejecutar la fase 3 (avance).
- [ ] Digitalizar la Fig. 2.7(e) y comparar los termopares.
- [ ] Averiguar la profundidad de los termopares en Stubblefield et al. (ref. [38] de la tesis).
- [ ] Pasar k del sustrato de constante (167 W/m·K) a la cúbica de la tesis (error < 10 % en el rango del proceso).
- [ ] Repetir con los casos *starved* (Vᵢ = 1.06 mm/s) y *overfed* (4.24 mm/s): cambiar `feedVelocity` en `toolProperties`.
- [ ] Probar `preparar_caso.sh` y `mallar.sh` en limpio: agrupan comandos probados uno a uno, pero no se han ejecutado como scripts.

## Licencia

GNU General Public License v3.0 (ver [`LICENSE`](LICENSE)). `afsdFoam` deriva de `chtMultiRegionFoam` de OpenFOAM (© OpenFOAM Foundation, © OpenCFD Ltd.), distribuido bajo GPL v3, así que este repositorio también lo es.

OpenFOAM® es una marca registrada de OpenCFD Ltd. Este proyecto no está aprobado ni respaldado por OpenCFD Ltd.
