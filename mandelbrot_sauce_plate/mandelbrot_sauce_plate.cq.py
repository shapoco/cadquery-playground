import cadquery as cq
import numpy as np

SCALE = 40
Z_SCALE = 20
BASE_R = 30
#PITCH_GOAL = 2
PITCH_GOAL = 0.4
BASE_R2 = BASE_R - PITCH_GOAL * 2
N_DIV = int(np.ceil(BASE_R2 * 2 / PITCH_GOAL))
PITCH = BASE_R2 * 2 / N_DIV
BASE_L = 3
ENG_DEPTH = 5
BASE_H = BASE_L + ENG_DEPTH
DEPTH = 5
THICKNESS = 5

AC = -0.5
BC = 0
AB_R = 1.5

A0 = AC - AB_R
A1 = AC + AB_R
B0 = BC - AB_R
B1 = BC + AB_R


def mandelbrot(a, b, max_iter):
    z = 0
    c = complex(a, b)
    for n in range(max_iter):
        if abs(z) > 2:
            return n
        z = z * z + c
    return max_iter


def build(ix0, ix1):
    ixc = (ix0 + ix1) // 2
    x0 = -BASE_R2 + BASE_R2 * 2 * ix0 / N_DIV
    x1 = -BASE_R2 + BASE_R2 * 2 * ix1 / N_DIV
    x_thick = x1 - x0
    x_middle = (x0 + x1) / 2
    a = AC + AB_R * x_middle / BASE_R2

    if ix0 + 1 == ix1:
        r = BASE_R2 - PITCH / 2
        iy0 = int(np.floor(-np.sqrt(r**2 - x_middle**2) / PITCH)) - 1
        iy1 = -iy0
        iy0 += int(N_DIV // 2)
        iy1 += int(N_DIV // 2)

        y0 = -BASE_R2 + BASE_R2 * 2 * iy0 / N_DIV
        y1 = -BASE_R2 + BASE_R2 * 2 * iy1 / N_DIV

        max_iter = 100
        z_last = -1
        verts = [(y0, 5)]
        for iy in range(iy0, iy1):
            y = -BASE_R2 + BASE_R2 * 2 * iy / N_DIV
            b = BC + AB_R * (y + PITCH / 2) / BASE_R2
            z = (mandelbrot(a, b, max_iter)) / max_iter
            z = 1 - (1 - z) ** 4
            z = -ENG_DEPTH * z
            if z != z_last:
                verts.append((y, z_last))
                verts.append((y, z))
            z_last = z
        verts.append((y1, z_last))
        verts.append((y1, 5))
        return (
            cq.Workplane("YZ")
            .polyline(verts)
            .close()
            .extrude(x_thick)
            .translate((x0, 0, 0))
        )
    else:
        sub0 = build(ix0, ixc)
        sub1 = build(ixc, ix1)
        return sub0.union(sub1)


dish_inner_verts = [
    (0, 0),
    (0, BASE_H),
    (BASE_R, BASE_H),
    (BASE_R + DEPTH * 1.5, BASE_H + DEPTH),
    (BASE_R + DEPTH * 1.5 + THICKNESS, BASE_H + DEPTH),
    (BASE_R + DEPTH * 0.5 - BASE_H, 0),
]

dish_inner = (
    cq.Workplane("YZ")
    .polyline(dish_inner_verts)
    .close()
    .revolve(360, (0, 0, 0), (0, 1, 0))
    .faces(">Z")
    .fillet(1)
    .faces("<Z")
    .fillet(5)
)

cutter = build(0, N_DIV).translate((0, 0, BASE_H + 0.4))

#show_object(cutter)

result = dish_inner.cut(cutter)

show_object(result)
result.export("mandelbrot.step")
