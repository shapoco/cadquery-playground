# https://htlab.net/computer/software/cadquery/cadquery-tutorial/

import cadquery as cq
import numpy as np

# 焦点距離 [mm]
F1 = 50
F2 = 200

# 円周方向の分割数
SPLIT = 12

# 反射面の半径
DISH_R = 75

# 反射面の厚み
DISH_T = 4

# 反射面の曲線の頂点数
N = (DISH_R + 1) // 2

# 中心の穴の半径
HOLE_R = 7.5

# 固定用の壁の高さ
WALL_H = 3

# 壁の厚み
WALL_T = 3

# マージン
MARGIN = 0.5

# 中心のディスク
DISC_R = 20
DISC_H = WALL_H * 2
DISC_HOLE_D = 3.5
DISC_SHOULDER_H = DISC_H - WALL_H

M3_CHAMFER = 3.5
M3_NUT_S = 5.4
M3_NUT_T = 3

# リングの太さ
RING_W = 3
RING_H = 5
RING_OFST_X = 10
RING_OFST_Y = 10

# 腕の数
NUM_ARMS = 6

# 腕の幅
ARM_W = 8

# center_x = 0
# center_y = (F2 + F1) / 2
radius_z = (F2 + F1) / 2
radius_y = np.sqrt(radius_z**2 - ((F2 - F1) / 2) ** 2)
angle_start = np.arcsin(HOLE_R / radius_y) * 180 / np.pi
angle_end = np.arcsin(DISH_R / radius_y) * 180 / np.pi

edge_norm_y = np.cos(angle_end * np.pi / 180)
edge_norm_z = -np.sin(angle_end * np.pi / 180) * radius_y / radius_z
edge_norm_d = np.sqrt(edge_norm_y**2 + edge_norm_z**2)
edge_norm_y /= edge_norm_d
edge_norm_z /= edge_norm_d

edge_x = np.sin(angle_end * np.pi / 180) * (radius_y + DISH_T)
edge_y = -np.cos(angle_end * np.pi / 180) * (radius_z + DISH_T)

hole_y = np.sin(angle_start * np.pi / 180) * radius_y
hole_z = -np.cos(angle_start * np.pi / 180) * radius_z

dish_offset = -hole_z + DISC_H + DISH_T

if False:

    # 楕円のパラメータ
    c = (F2 - F1) / 2
    a = (F2 + F1) / 2
    b = np.sqrt(a**2 - c**2)

    # 反射面の座標
    curve_x = np.linspace(HOLE_R, DISH_R, N)
    curve_y = []
    for x in curve_x:
        curve_y.append(-np.sqrt(a**2 * (1 - (x**2) / (b**2))))

    # 高さ調整
    curve_offset_y = DISC_H + DISH_T - min(curve_y)
    for i in range(N):
        curve_y[i] += curve_offset_y

    # 回転体のプロファイル
    dish_wall_verts = []
    edge_y = None
    edge_x = None
    for i in range(N):
        x = curve_x[i]
        y = curve_y[i]
        dish_wall_verts.append((x, y))

    # 厚みを持たせる
    for i in range(N - 1, -1, -1):
        x = curve_x[i]
        y = curve_y[i]
        norm_x = 0
        norm_y = 0
        if i == 0:
            # 穴の内側の壁は垂直に
            norm_x = 0
            norm_y = -1
        else:
            # 他の頂点は法線方向に
            if i > 0:
                norm_x -= curve_y[i - 1] - curve_y[i]
                norm_y += curve_x[i - 1] - curve_x[i]
            if i < N - 1:
                norm_x -= curve_y[i] - curve_y[i + 1]
                norm_y += curve_x[i] - curve_x[i + 1]
        d = np.sqrt(norm_y**2 + norm_x**2)
        norm_x = norm_x * DISH_T / d
        norm_y = norm_y * DISH_T / d
        x = x + norm_x
        y = y + norm_y
        dish_wall_verts.append((x, y))

        if i == N - 1:
            # 最も外側の座標
            edge_x = x
            edge_y = y

# 皿
dish_wall = (
    cq.Workplane("YZ")
    # .moveTo(HOLE_R, 0)
    .ellipseArc(
        radius_y,
        radius_z,
        angle_start - 90,
        angle_end - 90,
        sense=1,
        startAtCurrent=False,
    )
    .lineTo(edge_x, edge_y)
    .ellipseArc(
        radius_y + DISH_T,
        radius_z + DISH_T,
        angle_start - 90,
        angle_end - 90,
        sense=-1,
        startAtCurrent=True,
    )
    .close()
    .translate((0, 0, dish_offset))
)
show_object(dish_wall)
dish_wall = dish_wall.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0)).translate(
    (0, 0, dish_offset)
)

edge_y += dish_offset
# show_object(dish_wall)

# 穴のガイド
hole_guide_verts = [
    (HOLE_R, 0),
    (DISC_R + WALL_T, 0),
    (DISC_R + WALL_T, -DISH_T - WALL_H),
    (DISC_R, -DISH_T - WALL_H),
    (DISC_R, -DISH_T),
    (HOLE_R, -DISH_T),
]
hole_guide = (
    cq.Workplane("YZ")
    .polyline(hole_guide_verts)
    .close()
    .translate((0, 0, DISC_H + DISH_T))
    # .revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0))
)
show_object(hole_guide)
hole_guide = hole_guide.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0)).translate(
    (0, 0, DISC_H + DISH_T)
)

# リングのガイド
ring_guide_verts = [
    (0, 0),
    (-0.5, 0),
    (-RING_OFST_X, -RING_OFST_Y + 0.5),
    (-RING_OFST_X, -RING_OFST_Y),
    (-RING_W, -RING_OFST_Y),
    (-RING_W, -RING_OFST_Y + RING_H),
    (0, -RING_OFST_Y + RING_H),
]
for i in range(len(ring_guide_verts)):
    x, y = ring_guide_verts[i]
    x += edge_x
    y += edge_y
    ring_guide_verts[i] = (x, y)
ring_guide = cq.Workplane("YZ").polyline(ring_guide_verts).close()
show_object(ring_guide)
ring_guide = ring_guide.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0))

dish = dish_wall.union(hole_guide).union(ring_guide)

# 隣の皿とぶつからないように若干痩せさせる
dish = (
    dish.faces("<X")
    .workplane()
    .rect(edge_x * 2.5, edge_y * 2.5)
    .cutBlind(-MARGIN / 2)
    .faces(">X")
    .workplane()
    .rect(edge_x * 2.5, edge_y * 2.5)
    .cutBlind(-MARGIN / 2)
)

show_object(dish)

dish = dish.rotate((0, 0, 0), (0, 1, 0), 90).translate((0, 0, -MARGIN / 2))

dish.export("dish_single.step")

dishes = cq.Assembly()
for i in range(0, SPLIT):
    dishes.add(dish.translate((i * (6 + DISH_T), i * -(3 + DISH_T), 0)))
# show_object(dishes)

dishes.export("dishes.step")

disc_verts = [
    (0, 0),
    (DISC_R + WALL_T, 0),
    (DISC_R + WALL_T, DISC_SHOULDER_H),
    # (DISC_R - MARGIN / 2, DISC_SHOULDER_H),
    # (DISC_R - MARGIN / 2, DISC_SHOULDER_H + WALL_H),
    (DISC_R, DISC_SHOULDER_H),
    (DISC_R, DISC_SHOULDER_H + WALL_H),
    (M3_CHAMFER, DISC_SHOULDER_H + WALL_H),
    (0, DISC_SHOULDER_H + WALL_H - M3_CHAMFER),
]

disc = (
    cq.Workplane("YZ")
    .polyline(disc_verts)
    .close()
    .revolve(360, (0, 0, 0), (0, 1, 0))
    .faces(">Z")
    .workplane()
    .pushPoints([(0, 0, 0)])
    .circle(DISC_HOLE_D / 2)
    .cutBlind(-100)
    .faces("<Z")
    .workplane()
    # .polygon(6, M3_NUT_S * 2 / np.sqrt(3) + MARGIN / 2)
    .polygon(6, M3_NUT_S * 2 / np.sqrt(3))
    .cutBlind(-M3_NUT_T)
)

ring_outer_x = edge_x
ring_bottom_y = edge_y - RING_OFST_Y
ring_verts = [
    # (ring_outer_x - RING_W + MARGIN / 2, ring_bottom_y),
    (ring_outer_x - RING_W, ring_bottom_y),
    (ring_outer_x, ring_bottom_y),
    (ring_outer_x, ring_bottom_y + RING_H),
    # (ring_outer_x - RING_W + MARGIN / 2, ring_bottom_y + RING_H),
    (ring_outer_x - RING_W, ring_bottom_y + RING_H),
]

ring = (
    cq.Workplane("YZ").polyline(ring_verts).close().revolve(360, (0, 0, 0), (0, 1, 0))
)

ring_cutter = (
    ring.faces("<Z")
    .workplane()
    .rect(ARM_W + MARGIN, RING_W * 2)
    .extrude(-1 - MARGIN / 2, combine=False)
    .translate((0, ring_outer_x - RING_W / 2))
)
# show_object(ring_cutter)

for i in range(0, NUM_ARMS):
    ring = ring.cut(ring_cutter.rotate((0, 0, 0), (0, 0, 1), i * (360 / NUM_ARMS)))

show_object(ring)

ring = ring.rotate((0, 0, 0), (0, 1, 0), 180).translate((0, 0, ring_bottom_y + WALL_H))

# show_object(ring)

ring.export("ring.step")

p = WALL_T
arm_verts = [
    (DISC_R - WALL_T, 0),
    (ring_outer_x - ring_bottom_y + p, 0),
    (ring_outer_x, ring_bottom_y - p),
    (ring_outer_x, ring_bottom_y + 1),
    (ring_outer_x - RING_W, ring_bottom_y + 1),
    (ring_outer_x - RING_W, ring_bottom_y - 1),
    (ring_outer_x - ring_bottom_y - RING_W + WALL_T + 1, 0 + WALL_T),
    (DISC_R - WALL_T, WALL_T),
]

arm0 = (
    cq.Workplane("YZ")
    .polyline(arm_verts)
    .close()
    .extrude(ARM_W / 2, both=True)
    .faces("|Z")
    .faces("not (>Z or <Z)")
    .edges(">Y")
    .fillet(3)
    .faces("<Z")
    .edges(">Y")
    .fillet(5)
)

# プリント中に反らないように腕に溝を掘る
arm_cutter_verts = [
    (ring_outer_x - RING_W - 5, ring_bottom_y - 1),
    (ring_outer_x - ring_bottom_y - RING_W + WALL_T - 4, 0 + WALL_T),
    (ring_outer_x - ring_bottom_y - RING_W + WALL_T + 3, 0 + WALL_T),
    (ring_outer_x - RING_W, ring_bottom_y - 3),
    (ring_outer_x - RING_W, ring_bottom_y - 1),
]
arm_cutter = cq.Workplane("YZ").polyline(arm_cutter_verts).close().extrude(1, both=True)
arm0 = arm0.cut(arm_cutter)

frame = disc
for i in range(0, NUM_ARMS):
    frame = frame.union(arm0.rotate((0, 0, 0), (0, 0, 1), i * (360 / NUM_ARMS)))

# show_object(disc)
# show_object(ring)
# show_object(arm0)
# show_object(arm1)
# show_object(arm2)

show_object(frame)
frame.export("frame.step")
