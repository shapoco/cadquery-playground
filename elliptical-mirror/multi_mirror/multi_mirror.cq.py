import cadquery as cq

# import numpy as np
import math

# 反射板の数
NUM_REFLECTORS = 6

# 焦点間の距離 [mm]
ELLIPSE_FOCUS_DISTANCE = 40

# 楕円のパラメータ
ELLIPSE_Y_RADIUS = 40
ELLIPSE_X_RADIUS = math.sqrt((ELLIPSE_FOCUS_DISTANCE / 2) ** 2 + ELLIPSE_Y_RADIUS**2)
ELLIPSE_ANGLE = 45

# 嵌合用の隙間 [mm]
GAP = 0.2

# 固定用キャップの穴の間隔 [mm]
CAP_HOLE_DISTANCE = 10

# LED アームの取り付け用穴の間隔 [mm]
LED_ARM_HOLE_DISTANCE = 15

STEP_OUT_DIR = "./step"

class Reflector:
    def __init__(self):
        inner_radius = 10
        outer_radius = 55
        height = 65

        inner_sin = inner_radius * math.sin(math.pi / NUM_REFLECTORS)
        inner_cos = inner_radius * math.cos(math.pi / NUM_REFLECTORS)
        solid = (
            cq.Workplane("XY")
            .ellipseArc(
                outer_radius,
                outer_radius,
                -180 / NUM_REFLECTORS,
                180 / NUM_REFLECTORS,
                startAtCurrent=False,
            )
            .lineTo(inner_cos, inner_sin)
            .lineTo(inner_cos, -inner_sin)
            .close()
            .extrude(height)
        )

        # 楕円面を彫り込む
        cutter = (
            cq.Workplane("XY")
            .box(400, 400, 200, centered=(True, True, False))
            .translate((0, 0, -200))
        )
        cutter = cutter.union(
            cq.Workplane("XY")
            .ellipseArc(
                ELLIPSE_X_RADIUS,
                ELLIPSE_Y_RADIUS,
                0,
                180,
                sense=-1,
                startAtCurrent=False,
            )
            .close()
            .revolve(360, (0, 0, 0), (1, 0, 0))
            .translate((ELLIPSE_FOCUS_DISTANCE / 2, 0, 0))
        )
        cutter = cutter.rotate((0, 1, 0), (0, 0, 0), ELLIPSE_ANGLE)
        solid = solid.cut(cutter)

        # LED 取り付け用の溝を掘る
        x = outer_radius - 5 - GAP
        y = 15 + GAP
        verts = [
            (x, -y),
            (x, y),
            (x + 10, y + 10),
            (x + 10, -y - 10),
        ]
        cutter = (
            cq.Workplane("XY")
            .polyline(verts)
            .close()
            .extrude(height + 10)
        )
        solid = solid.cut(cutter)

        # 天井を斜めにする
        verts = [
            (outer_radius - 10, height + 10),
            (outer_radius, height + 10),
            (outer_radius, height - 10),
            (outer_radius - 10, height),
        ]
        cutter = (
            cq.Workplane("XZ").polyline(verts).close().extrude(outer_radius, both=True)
        )
        solid = solid.cut(cutter)

        # 穴開け
        solid = (
            solid.faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (25, -CAP_HOLE_DISTANCE / 2),
                    (25, CAP_HOLE_DISTANCE / 2),
                ]
            )
            .circle(1)
            .cutBlind(-5)
            .faces(">>X[-2]")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (-LED_ARM_HOLE_DISTANCE / 2, height - 10),
                    (LED_ARM_HOLE_DISTANCE / 2, height - 10),
                ]
            )
            .circle(1)
            .cutBlind(-5)
        )

        # 面取り
        solid = (
           solid.faces("<X or <<Y[-2] or >>Y[-2] or <<Z[-2]")
           .chamfer(1)
           .edges("%circle")
           .chamfer(0.5)
        )

        self.solid = solid

ref = Reflector()

show_object(ref.solid)

# 焦点の位置 (参考用)
show_object(
    cq.Workplane("XY")
    .box(3, 3, 3)
    .translate((ELLIPSE_FOCUS_DISTANCE, 0, 0))
    .rotate((0, 1, 0), (0, 0, 0), ELLIPSE_ANGLE)
)
show_object(
    cq.Workplane("XY")
    .box(3, 3, 3)
    .translate((0, 0, 0))
    .rotate((0, 1, 0), (0, 0, 0), ELLIPSE_ANGLE)
)

ref.solid.export(f"{STEP_OUT_DIR}/trial.step")
