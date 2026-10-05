"""
Parámetros de la geometría y la malla AFSD (Kincaid 2022, cap. 2).
Todas las longitudes en mm; la malla se escala a metros al exportar.

Marco de referencia: herramienta fija en el origen, avanza en +x.
z = 0 es la cara superior del sustrato.
"""
import math

# --- Geometría (tesis) -------------------------------------------------------
L_DOMAIN = 152.4     # longitud del sustrato / dominio en x            (cap. 2)
W_DOMAIN = 76.2      # anchura en y                                     (cap. 2)
T_SUB    = 6.35      # espesor del sustrato                             (cap. 2)
T_ANVIL  = 12.7      # espesor del yunque                               (cap. 4)
H_LAYER  = 1.5       # altura de la herramienta sobre el sustrato       (cap. 2)
H_BUILD  = 50.0      # altura total de la región build                  (cap. 4)
D_TOOL   = 38.1      # diámetro de la herramienta                       (cap. 2)
ROD_SIDE = 9.35      # lado de la barra cuadrada                        (cap. 3)

R_TOOL  = D_TOOL / 2
R_INLET = math.sqrt(ROD_SIDE**2 / math.pi)      # círculo de igual área (cap. 2)

# O-grid bajo la herramienta: cuadrado central dentro del círculo de entrada.
# Radio de las esquinas del cuadrado como fracción de R_INLET.
SQUARE_CORNER_FRAC = 0.55
R_SQUARE_CORNER = SQUARE_CORNER_FRAC * R_INLET

# --- Malla -------------------------------------------------------------------
N_CIRC      = 12     # divisiones por cuarto de circunferencia (y por lado del cuadrado)
N_RAD_INLET = 5      # radiales entre el cuadrado y el borde de la entrada
N_RAD_TOOL  = 14     # radiales entre la entrada y el borde de la herramienta (cap. 5)
N_GAP       = 9      # celdas en z entre herramienta y sustrato           (cap. 5)

# Capas graduadas: primera celda junto a la interfaz más cercana
N_UP,    DZ1_UP    = 20, 0.25   # build por encima de la herramienta
N_SUB,   DZ1_SUB   = 10, 0.25   # sustrato (fino junto a z = 0)
N_ANVIL, DZ1_ANVIL = 6,  1.0    # yunque   (fino junto al sustrato)

# Tamaño de celda en planta fuera de la herramienta
SIZE_PATH = 2.0      # banda de la trayectoria |y| < Y_PATH
Y_PATH    = 25.0
SIZE_FAR  = 6.0      # resto
RECOMBINE_OUTER = False   # False: triángulos -> prismas (como la tesis)

MESH_FILE  = "afsd.msh"
BREP_FILE  = "planform.brep"
STEP_3D    = "domain_3d.step"   # solo para visualizar
