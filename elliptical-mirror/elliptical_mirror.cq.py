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

# ノッチの高さ
NOTCH_H = 2

# 中心のディスクの寸法
DISC_R = 20
DISC_H = WALL_H * 2 + NOTCH_H
DISC_HOLE_D = 3.5
DISC_SHOULDER_H = DISC_H - WALL_H

# リングの寸法
RING_W = 3
RING_H = 5
RING_OFST_X = 10
RING_OFST_Y = 10

# 腕の寸法
NUM_ARMS = 6
ARM_W = 8

# M3ネジ寸法
M3_CHAMFER = 3.5
M3_NUT_S = 5.4
M3_NUT_T = 3

OUT_DIR = "./step"

# 警告がウザいので自前の wrapper をかます
def show_obj(obj):
    show_object(obj)


# お皿の部分
class Dish:
    def __init__(self):
        self.solid = None

        # 楕円のパラメータ
        radius_z = (F2 + F1) / 2
        radius_y = np.sqrt(radius_z**2 - ((F2 - F1) / 2) ** 2)
        angle_start = np.arcsin(HOLE_R / radius_y) * 180 / np.pi
        angle_end = np.arcsin(DISH_R / radius_y) * 180 / np.pi

        # 穴の座標
        hole_y = np.sin(angle_start * np.pi / 180) * radius_y
        hole_z = -np.cos(angle_start * np.pi / 180) * radius_z

        # 皿の外縁部の座標
        edge_norm_x = np.cos(angle_end * np.pi / 180)
        edge_norm_y = -np.sin(angle_end * np.pi / 180) * radius_y / radius_z
        d = np.sqrt(edge_norm_x**2 + edge_norm_y**2)
        edge_norm_x /= d
        edge_norm_y /= d
        edge_y = np.sin(angle_end * np.pi / 180) * (radius_y + DISH_T)
        edge_z = -np.cos(angle_end * np.pi / 180) * (radius_z + DISH_T)

        # 高さ調整
        z_offset = -hole_z + DISC_H + DISH_T

        # 皿
        wall_profile = (
            cq.Workplane("YZ")
            .ellipseArc(
                radius_y,
                radius_z,
                angle_start - 90,
                angle_end - 90,
                sense=1,
                startAtCurrent=False,
            )
            .lineTo(edge_y, edge_z)
            .ellipseArc(
                radius_y + DISH_T,
                radius_z + DISH_T,
                angle_start - 90,
                angle_end - 90,
                sense=-1,
                startAtCurrent=True,
            )
            .close()
            .translate((0, 0, z_offset))
        )
        wall = wall_profile.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0)).translate(
            (0, 0, z_offset)
        )

        edge_z += z_offset

        # 穴のガイド
        z = DISC_H + DISH_T
        inner_guide_verts = [
            (HOLE_R, z),
            (DISC_R + WALL_T, z),
            (DISC_R + WALL_T, z - DISH_T - WALL_H),
            (DISC_R, z - DISH_T - WALL_H),
            (DISC_R, z - DISH_T),
            (HOLE_R, z - DISH_T),
        ]
        inner_guide_profile = cq.Workplane("YZ").polyline(inner_guide_verts).close()
        hole_guide = inner_guide_profile.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0))

        # 治具作成のために頂点抽出
        jig_verts0 = [
            inner_guide_verts[1],
            inner_guide_verts[2],
            inner_guide_verts[3],
        ]

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
            x += edge_y
            y += edge_z
            ring_guide_verts[i] = (x, y)
        outer_guide_profile = cq.Workplane("YZ").polyline(ring_guide_verts).close()
        ring_guide = outer_guide_profile.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0))

        # 治具作成のために頂点抽出
        jig_verts1 = [
            ring_guide_verts[6],
            ring_guide_verts[5],
            ring_guide_verts[4],
            ring_guide_verts[3],
        ]

        model_without_margin = wall.union(hole_guide).union(ring_guide)

        # 隣の皿とぶつからないように若干痩せさせる
        self.solid = (
            model_without_margin.faces("<X")
            .workplane()
            .rect(edge_y * 2.5, edge_z * 2.5)
            .cutBlind(-MARGIN / 2)
            .faces(">X")
            .workplane()
            .rect(edge_y * 2.5, edge_z * 2.5)
            .cutBlind(-MARGIN / 2)
        )

        self.model_without_margin = model_without_margin
        self.edge_x = edge_y
        self.edge_y = edge_z
        self.jig_verts0 = jig_verts0
        self.jig_verts1 = jig_verts1
        self.wall_profile = wall_profile
        self.inner_guide_profile = inner_guide_profile
        self.outer_guide_profile = outer_guide_profile

    def export(self):
        # 印刷用にビルドボードに立てる
        single = self.solid.rotate((0, 0, 0), (0, 1, 0), 90).translate(
            (0, 0, -MARGIN / 2)
        )
        single.export(f"{OUT_DIR}/dish_single.step")

        # まとめて印刷する用に並べる
        array = cq.Assembly()
        for i in range(0, SPLIT):
            array.add(single.translate((i * (6 + DISH_T), i * -(3 + DISH_T), 0)))
        array.export(f"{OUT_DIR}/dish_array.step")


class DishJig:
    def __init__(self, dish: Dish):

        v00 = dish.jig_verts0[0]
        v01 = dish.jig_verts0[1]
        v02 = dish.jig_verts0[2]
        v03 = (v02[0], v02[1] - WALL_T)
        v06 = dish.jig_verts1[0]
        v07 = dish.jig_verts1[1]
        v08 = dish.jig_verts1[2]
        v09 = dish.jig_verts1[3]

        # 接地面の向きと角度
        dy = v06[0] - v03[0]
        dz = v06[1] - v03[1]
        angle = np.arctan2(dz, dy) * 180 / np.pi

        # 接地面の法線
        norm_y = dz
        norm_z = -dy
        norm_d = np.sqrt(norm_y**2 + norm_z**2)
        norm_y /= norm_d
        norm_z /= norm_d

        # 接地面に厚みを持たせる
        h = 10
        v04 = (v03[0] + norm_y * h, v03[1] + norm_z * h)
        v05 = (v06[0] + norm_y * h, v06[1] + norm_z * h)

        # 皿の形状に合わせて中点を作成
        h = WALL_T
        v10 = ((v09[0] + v00[0]) / 2 + norm_y * h, (v09[1] + v00[1]) / 2 + norm_z * h)

        # 回転体のプロファイル
        verts = [v00, v01, v02, v03, v04, v05, v06, v07, v08, v09, v10]
        profile = cq.Workplane("YZ").polyline(verts).close()
        solid = (
            profile.revolve(360 / SPLIT, (0, 0, 0), (0, 1, 0))
            .faces(">X")
            .extrude(1.5)
            .rotate((0, 0, 0), (0, 0, 1), -360 / SPLIT)
            .faces("<X")
            .extrude(-1.5)
            .rotate((0, 0, 0), (0, 0, 1), 360 / SPLIT)
            .cut(dish.model_without_margin)
        )

        # プラットフォームと平行にする
        r_axis0 = (0, v03[0], v03[1])
        r_axis1 = (
            r_axis0[1] * np.sin(-2 * np.pi / SPLIT),
            r_axis0[1] * np.cos(-2 * np.pi / SPLIT),
            r_axis0[2],
        )
        solid = solid.rotate(r_axis0, r_axis1, angle).translate(
            (0, 0, -v03[1] + WALL_T * 2)
        )

        # 原点付近へ移動
        solid = solid.rotate((0, 0, 0), (0, 0, 1), -180 / SPLIT).translate(
            (0, -DISH_R * 3 / 4, 0)
        )

        # プラットフォームより下を除去して平らにする
        floor_cutter = (
            cq.Workplane("XY").rect(-DISH_R * 2, DISH_R * 2).extrude(-WALL_T * 10)
        )
        solid = solid.cut(floor_cutter)

        hole_verts = [
            (-12, DISH_R * 0.3),
            (12, DISH_R * 0.3),
            (4, -DISH_R * 0.3),
            (-4, -DISH_R * 0.3),
        ]
        hole_cutter = (
            cq.Workplane("XY")
            .polyline(hole_verts)
            .close()
            .extrude(50)
            .edges("|Z")
            .fillet(3)
            .translate((0, 0, 3))
        )
        solid = (
            solid.cut(hole_cutter)
            # ネジ穴
            .faces("<Z")
            .workplane()
            .pushPoints([(0, -15), (0, 15)])
            .circle(3.5 / 2)
            .cutBlind(-10)
            # 皿ネジ用の面取り
            .faces("+Z")
            .edges()[-2:]
            .chamfer(2)
            # 何かに使うかもしれない中央の穴
            .faces("<Z")
            .workplane()
            .pushPoints([(0, 0)])
            .circle(5)
            .cutBlind(-10)
        )

        self.profile = profile
        self.solid = solid

    def export(self):
        self.solid.export(f"{OUT_DIR}/dish_jig.step")


class Ring:
    def __init__(self, dish: Dish):
        # リングの原型
        outer_x = dish.edge_x
        bottom_y = dish.edge_y - RING_OFST_Y
        verts = [
            (outer_x - RING_W, bottom_y),
            (outer_x, bottom_y),
            (outer_x, bottom_y + RING_H),
            (outer_x - RING_W, bottom_y + RING_H),
        ]
        solid = (
            cq.Workplane("YZ")
            .polyline(verts)
            .close()
            .revolve(360, (0, 0, 0), (0, 1, 0))
        )

        # 腕を取り付ける溝を彫る
        cutter = (
            solid.faces("<Z")
            .workplane()
            .rect(ARM_W + MARGIN, RING_W * 2)
            .extrude(-1 - MARGIN / 2, combine=False)
            .translate((0, outer_x - RING_W / 2))
        )
        for i in range(0, NUM_ARMS):
            solid = solid.cut(cutter.rotate((0, 0, 0), (0, 0, 1), i * (360 / NUM_ARMS)))

        self.solid = solid
        self.outer_x = outer_x
        self.bottom_y = bottom_y

    def export(self):
        flipped = self.solid.rotate((0, 0, 0), (0, 1, 0), 180).translate(
            (0, 0, self.bottom_y + RING_H)
        )
        flipped.export(f"{OUT_DIR}/ring.step")


class Frame:
    def __init__(self, ring: Ring):

        # 中央の円盤
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

        # 腕
        p = WALL_T
        arm_verts = [
            (DISC_R - WALL_T, 0),
            (ring.outer_x - ring.bottom_y + p, 0),
            (ring.outer_x, ring.bottom_y - p),
            (ring.outer_x, ring.bottom_y + 1),
            (ring.outer_x - RING_W, ring.bottom_y + 1),
            (ring.outer_x - RING_W, ring.bottom_y - 1),
            (ring.outer_x - ring.bottom_y - RING_W + WALL_T + 1, 0 + WALL_T),
            (DISC_R - WALL_T, WALL_T),
        ]

        # 腕のテンプレート
        arm = cq.Workplane("YZ").polyline(arm_verts).close()
        arm = arm.extrude(ARM_W / 2, both=True)

        # 腕がプリント中に反らないためのノッチ (意味無いかも)
        q = ring.bottom_y * 1 / 3
        notch_verts = [
            (DISC_R - WALL_T, WALL_T + NOTCH_H),
            (DISC_R + WALL_T + 2, WALL_T + NOTCH_H),
            (ring.outer_x - RING_W - NOTCH_H * 3 / 2 - q, ring.bottom_y - 1 - q),
            (ring.outer_x - RING_W - NOTCH_H * 3 / 2, ring.bottom_y - 1),
            (ring.outer_x - RING_W, ring.bottom_y - 1),
            (ring.outer_x - RING_W, ring.bottom_y - 2),
            (ring.outer_x - ring.bottom_y + 1, WALL_T - 1),
            (DISC_R - WALL_T, WALL_T - 1),
        ]
        notch = cq.Workplane("YZ").polyline(notch_verts).close()
        notch = notch.extrude(1, both=True)
        arm = arm.union(notch)

        # 円盤に腕をくっつける
        solid = disc
        for i in range(0, NUM_ARMS):
            solid = solid.union(arm.rotate((0, 0, 0), (0, 0, 1), i * (360 / NUM_ARMS)))

        self.solid = solid

    def export(self):
        self.solid.export(f"{OUT_DIR}/frame.step")


dish = Dish()
jig = DishJig(dish)
ring = Ring(dish)
frame = Frame(ring)

show_object(dish.wall_profile)
show_object(dish.inner_guide_profile)
show_object(dish.outer_guide_profile)
show_object(dish.solid)
show_object(jig.profile)
show_object(jig.solid.translate((DISH_R, 0, 0)))
show_object(ring.solid)
show_object(frame.solid)

dish.export()
jig.export()
ring.export()
frame.export()
