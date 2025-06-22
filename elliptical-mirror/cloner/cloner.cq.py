import cadquery as cq
import math

# 焦点間の距離 [mm]
F_DISTANCE = 100

# 鏡の設置高さ [mm]
MIRROR_HEIGHT = 150

# 鏡の一辺の長さ [mm]
MIRROR_WIDTH = 80

# 鏡の厚さ
REF_T = 5

# 楕円のパラメータ
INNER_RY = MIRROR_HEIGHT
INNER_RX = math.sqrt((F_DISTANCE / 2) ** 2 + INNER_RY**2)
print(f"RADIUS_X: {INNER_RX}")
OUTER_RY = INNER_RY + REF_T
OUTER_RX = INNER_RX + REF_T

# 支柱 (針金) の太さ
WIRE_DIAMETER = 2

# 出力先ディレクトリ
STEP_OUT_DIR = "./step"
STL_OUT_DIR = "./stl"


class Reflector:
    def __init__(self):
        # 半天球を覆う鏡
        solid = (
            cq.Workplane("XY")
            .ellipseArc(INNER_RX, INNER_RY, 0, 180, sense=-1, startAtCurrent=False)
            .lineTo(OUTER_RX, 0)
            .ellipseArc(OUTER_RX, OUTER_RY, 0, 180, sense=1, startAtCurrent=True)
            .close()
            .revolve(180, (0, 0, 0), (1, 0, 0))
        )

        # 鏡の大きさの箱で切り取る
        intersector = (
            cq.Workplane("XY")
            .box(MIRROR_WIDTH, MIRROR_WIDTH, 100)
            .translate((0, 0, MIRROR_HEIGHT))
        )
        solid = solid.intersect(intersector)

        # 首の部分
        neck = (
            cq.Workplane("XY")
            .box(MIRROR_WIDTH * 2 / 3, MIRROR_WIDTH / 2, 20)
            .translate((0, -MIRROR_WIDTH / 4, MIRROR_HEIGHT + 3 - 5))
        )
        # 鏡の前面にはみ出さないように切り取る
        body_cutter = (
            cq.Workplane("XY")
            .ellipseArc(
                OUTER_RX - 1, OUTER_RY - 1, 0, 180, sense=-1, startAtCurrent=False
            )
            .close()
            .revolve(180, (0, 0, 0), (1, 0, 0))
        )
        neck = neck.cut(body_cutter)
        solid = solid.union(neck)

        # 原点付近に移動
        solid = solid.translate((0, 0, -MIRROR_HEIGHT - 1))

        solid = (
            # 支柱用の穴開け
            solid.faces("<Y")
            .workplane()
            .pushPoints(
                [
                    (-15, 0),
                    (15, 0),
                ]
            )
            .circle(WIRE_DIAMETER / 2 + 0.25)
            .cutBlind(-20)
            # 丸め
            .edges("|Z")
            .fillet(3)
            .edges(">Z")
            .edges("not <Y")
            .fillet(3)
        )

        # 3Dプリント向けに垂直に立てる
        solid = solid.rotate(
            (0, -MIRROR_WIDTH / 2, 0), (1, -MIRROR_WIDTH / 2, 0), 90
        ).translate((0, MIRROR_WIDTH / 2, 0))

        self.solid = solid


class Base:
    def __init__(self):
        solid = (
            cq.Workplane("XY")
            # 床部分
            .box(MIRROR_WIDTH, MIRROR_WIDTH * 2 / 3, 5, centered=(True, False, False))
            .faces(">Z")
            .workplane()
            .polyline(
                [
                    (-MIRROR_WIDTH / 2 + 15, 15),
                    (-MIRROR_WIDTH / 2 + 15, MIRROR_WIDTH),
                    (MIRROR_WIDTH / 2 - 15, MIRROR_WIDTH),
                    (MIRROR_WIDTH / 2 - 15, 15),
                ]
            )
            .close()
            .cutBlind(-100)
            
            # 支柱取り付け部分の押し出し
            .faces(">Z")
            .workplane()
            .polyline(
                [
                    (-25, 0),
                    (-25, 10),
                    (25, 10),
                    (25, 0),
                ]
            )
            .close()
            .extrude(20)
            
            # 丸め
            .edges("|Z")
            .fillet(3)
            .edges(">Z")
            .fillet(3)
            
            # 支柱用の穴開け
            .faces(">Z")
            .workplane()
            .pushPoints(
                [
                    (-15, 5),
                    (15, 5),
                ]
            )
            .circle(WIRE_DIAMETER / 2 + 0.25)
            .cutBlind(-20)
            
        )
        self.solid = solid


ref = Reflector()
base = Base()

show_object(ref.solid, name="Reflector")
show_object(base.solid.translate((0, 20, 0)), name="Base")

base.solid.export(f"{STEP_OUT_DIR}/base.step")
base.solid.export(f"{STL_OUT_DIR}/base.stl")
