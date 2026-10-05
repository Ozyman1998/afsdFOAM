"""
Geometría paramétrica con CadQuery.

Genera la planta (z = 0) particionada para el mallado por extrusión:
  - cuadrado central + 4 sectores hasta el círculo de entrada  (O-grid)
  - 4 sectores entre la entrada y el borde de la herramienta    (O-grid)
  - rectángulo exterior menos el disco de la herramienta        (no estructurado)

Salidas:
  planform.brep    -> lo lee mesh.py
  domain_3d.step   -> sólidos build / substrate / anvil, solo para visualizar
"""
import math

import cadquery as cq
from cadquery import Edge, Face, Vector, Wire

from params import *


def polar(r, ang_deg, z=0.0):
    a = math.radians(ang_deg)
    return Vector(r * math.cos(a), r * math.sin(a), z)


def sector(r_in_pts, r_out, a0, a1):
    """Cara entre dos esquinas interiores (recta o arco) y un arco exterior."""
    p0, p1, inner_edge = r_in_pts
    q0, q1 = polar(r_out, a0), polar(r_out, a1)
    qm = polar(r_out, 0.5 * (a0 + a1))
    edges = [
        Edge.makeLine(p0, q0),
        Edge.makeThreePointArc(q0, qm, q1),
        Edge.makeLine(q1, p1),
        inner_edge,
    ]
    return Face.makeFromWires(Wire.assembleEdges(edges))


def planform_faces():
    faces = []
    corners = [45, 135, 225, 315, 405]

    # Cuadrado central
    sq = [polar(R_SQUARE_CORNER, a) for a in corners[:4]]
    faces.append(Face.makeFromWires(Wire.makePolygon(sq, close=True)))

    for k in range(4):
        a0, a1 = corners[k], corners[k + 1]

        # Sector cuadrado -> círculo de entrada
        p0, p1 = polar(R_SQUARE_CORNER, a0), polar(R_SQUARE_CORNER, a1)
        faces.append(sector((p0, p1, Edge.makeLine(p1, p0)), R_INLET, a0, a1))

        # Sector entrada -> borde de herramienta
        p0, p1 = polar(R_INLET, a0), polar(R_INLET, a1)
        pm = polar(R_INLET, 0.5 * (a0 + a1))
        faces.append(
            sector((p0, p1, Edge.makeThreePointArc(p1, pm, p0)), R_TOOL, a0, a1)
        )

    # Exterior: rectángulo menos el disco de la herramienta
    rect = Wire.makePolygon(
        [
            Vector(-L_DOMAIN / 2, -W_DOMAIN / 2, 0),
            Vector(L_DOMAIN / 2, -W_DOMAIN / 2, 0),
            Vector(L_DOMAIN / 2, W_DOMAIN / 2, 0),
            Vector(-L_DOMAIN / 2, W_DOMAIN / 2, 0),
        ],
        close=True,
    )
    # Agujero hecho con los mismos 4 arcos que los sectores: un círculo
    # completo tendría su vértice de costura en 0° y partiría un sector.
    hole = Wire.assembleEdges(
        [
            Edge.makeThreePointArc(
                polar(R_TOOL, a), polar(R_TOOL, a + 45), polar(R_TOOL, a + 90)
            )
            for a in corners[:4]
        ]
    )
    faces.append(Face.makeFromWires(rect, [hole]))
    return faces


def domain_3d():
    """Sólidos de referencia (no se mallan): útiles para revisar la geometría."""
    box = lambda z0, h: (
        cq.Workplane("XY").workplane(offset=z0).rect(L_DOMAIN, W_DOMAIN).extrude(h)
    )
    tool = (
        cq.Workplane("XY").workplane(offset=H_LAYER)
        .circle(R_TOOL).extrude(H_BUILD - H_LAYER)
    )
    build = box(0, H_BUILD).cut(tool)
    substrate = box(-T_SUB, T_SUB)
    anvil = box(-T_SUB - T_ANVIL, T_ANVIL)
    return (
        cq.Assembly()
        .add(build, name="build", color=cq.Color(0.6, 0.8, 1.0, 0.3))
        .add(substrate, name="substrate", color=cq.Color(0.2, 0.4, 0.9))
        .add(anvil, name="anvil", color=cq.Color(0.8, 0.2, 0.2))
        .add(tool, name="tool", color=cq.Color(0.2, 0.7, 0.2))
    )


if __name__ == "__main__":
    faces = planform_faces()
    cq.Compound.makeCompound(faces).exportBrep(BREP_FILE)
    domain_3d().export(STEP_3D)

    area = sum(f.Area() for f in faces)
    print(f"R_inlet = {R_INLET:.3f} mm   ({len(faces)} caras en planta)")
    print(f"Área planta = {area:.2f} mm²  (esperado {L_DOMAIN * W_DOMAIN:.2f})")
    print(f"Escrito: {BREP_FILE}, {STEP_3D}")
