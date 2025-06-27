import cadquery as cq

# import numpy as np
import math

# 反射板の数
NUM_REFLECTORS = 6

# 反射板の寸法 [mm]
REF_INNER_RADIUS = 10
REF_OUTER_RADIUS = 60
REF_HEIGHT = 65

# 焦点間の距離 [mm]
ELLIPSE_FOCUS_DISTANCE = 40

# 楕円のパラメータ
ELLIPSE_Y_RADIUS = 40
ELLIPSE_X_RADIUS = math.sqrt((ELLIPSE_FOCUS_DISTANCE / 2) ** 2 + ELLIPSE_Y_RADIUS**2)
ELLIPSE_ANGLE = 45

# 嵌合用の隙間 [mm]
GAP = 0.2

# 面取りのサイズ [mm]
CHAMFER = 0.5

# 固定用キャップの穴の間隔 [mm]
REF_INNER_MOUNT_X = 15
REF_OUTER_MOUNT_X = 30
REF_OUTER_MOUNT_Y = 10

# 固定用キャップの直径 [mm]
CAP_DIAMETER = (REF_OUTER_MOUNT_X + 12.5) * 2

# 固定用キャップの穴の間隔 [mm]
CAP_MOUNT_X = 15
CAP_MOUNT_Y = 20

# LED アームの取り付け用穴の間隔 [mm]
LED_ARM_HOLE_DISTANCE = 10

ARM_HOLDER_HEIGHT = 30

STEP_OUT_DIR = "./step"


class Reflector:
    def __init__(self):
        inner_sin = REF_INNER_RADIUS * math.sin(math.pi / NUM_REFLECTORS)
        inner_cos = REF_INNER_RADIUS * math.cos(math.pi / NUM_REFLECTORS)
        solid = (
            cq.Workplane("XY")
            .ellipseArc(
                REF_OUTER_RADIUS,
                REF_OUTER_RADIUS,
                -180 / NUM_REFLECTORS,
                180 / NUM_REFLECTORS,
                startAtCurrent=False,
            )
            .lineTo(inner_cos, inner_sin)
            .lineTo(inner_cos, -inner_sin)
            .close()
            .extrude(REF_HEIGHT)
        )

        # 楕円面を彫り込む
        cutter = (
            cq.Workplane("XY")
            .box(400, 400, 200, centered=(True, True, False))
            .translate((0, 0, -200 - 10))
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
        x = REF_OUTER_RADIUS - 5
        y = 10
        verts = [
            (x, -y),
            (x, y),
            (x + 10, y + 10),
            (x + 10, -y - 10),
        ]
        cutter = cq.Workplane("XY").polyline(verts).close().extrude(REF_HEIGHT + 10)
        solid = solid.cut(cutter)

        # 天井を斜めにする
        verts = [
            (REF_OUTER_RADIUS - 10, REF_HEIGHT + 10),
            (REF_OUTER_RADIUS, REF_HEIGHT + 10),
            (REF_OUTER_RADIUS, REF_HEIGHT - 10),
            (REF_OUTER_RADIUS - 10, REF_HEIGHT),
        ]
        cutter = (
            cq.Workplane("XZ")
            .polyline(verts)
            .close()
            .extrude(REF_OUTER_RADIUS, both=True)
        )
        solid = solid.cut(cutter)

        # 穴開け
        solid = (
            solid.faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (REF_INNER_MOUNT_X, 0),
                    (REF_OUTER_MOUNT_X, -REF_OUTER_MOUNT_Y),
                    (REF_OUTER_MOUNT_X, REF_OUTER_MOUNT_Y),
                ]
            )
            .circle(2.5 / 2)
            .cutBlind(-5)
            .faces(">>X[-3]")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (0, REF_HEIGHT - 10),
                    (0, REF_HEIGHT - 10 - LED_ARM_HOLE_DISTANCE),
                ]
            )
            .circle(2.5 / 2)
            .cutBlind(-5)
        )

        # 面取り
        solid = (
            solid
            # .faces("<X or <<Y[-2] or >>Y[-2] or <<Z[-2]")
            .faces("<X or <<Y[-2] or >>Y[-2] or <<Z[-3]")
            .chamfer(CHAMFER)
            .edges("%circle")
            .chamfer(CHAMFER)
        )

        self.solid = solid


class Cap:
    def __init__(self):
        solid = (
            cq.Workplane("XY")
            .polygon(NUM_REFLECTORS, CAP_DIAMETER)
            .extrude(5)
            .faces(">Z")
            .workplane()
            .polygon(NUM_REFLECTORS, REF_INNER_RADIUS * 2)
            .cutBlind(-50)
            .rotate((0, 0, 0), (0, 0, 1), 180 / NUM_REFLECTORS)
            .translate((0, 0, REF_HEIGHT))
        )

        # 穴開け
        hole_points_template = [
            (REF_INNER_MOUNT_X, 0),
            (REF_OUTER_MOUNT_X, -REF_OUTER_MOUNT_Y),
            (REF_OUTER_MOUNT_X, REF_OUTER_MOUNT_Y),
        ]
        for i in range(NUM_REFLECTORS):
            rad = math.radians(i * (360 / NUM_REFLECTORS))
            sin = math.sin(rad)
            cos = math.cos(rad)

            hole_points = []
            for j in range(len(hole_points_template)):
                x, y = hole_points_template[j]
                hole_points.append((x * cos - y * sin, x * sin + y * cos))

            solid = (
                solid.faces(">Z")
                .workplane(origin=(0, 0, 0))
                .pushPoints(hole_points)
                .circle(3.5 / 2)
                .cutBlind(-30)
            )

        # 面取り
        solid = (
            solid.edges("not %circle")
            .chamfer(CHAMFER)
            .edges(">Z and %circle")
            .chamfer(2)
            .edges("<Z and %circle")
            .chamfer(CHAMFER)
        )

        # マウント用の穴開け
        hole_points = [
            (CAP_MOUNT_X, -CAP_MOUNT_Y),
            (CAP_MOUNT_X, CAP_MOUNT_Y),
            (-CAP_MOUNT_X, -CAP_MOUNT_Y),
            (-CAP_MOUNT_X, CAP_MOUNT_Y),
        ]
        solid = (
            solid.faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints(hole_points)
            .circle(3.5 / 2)
            .cutBlind(-30)
            # 面取り
            .edges(">Z and %circle")
            .edges(">>Y[-3] or <<Y[-3]")
            .chamfer(CHAMFER)
            .edges("<Z and %circle")
            .edges(">>Y[-3] or <<Y[-3]")
            .chamfer(2)
        )

        self.solid = solid


class ArmHolder:
    def __init__(self):
        verts = [
            (0, -10),
            (2.5, -10),
            (2.5, -7.5),
            (5, -7.5),
            (5, 7.5),
            (2.5, 7.5),
            (2.5, 10),
            (0, 10),
        ]
        solid = (
            cq.Workplane("XY")
            .polyline(verts)
            .close()
            .extrude(ARM_HOLDER_HEIGHT)
            # 穴開け
            .faces(">X")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 5), (0, 15), (0, 25)])
            .circle(3.5 / 2)
            .cutBlind(-50)
            # 面取り
            .edges("%circle")
            .edges("(<X and (>Z or >>Z[-2])) or (>X and <Z)")
            .chamfer(CHAMFER)
            .edges("%circle")
            .edges(">X and (>Z or >>Z[-2]) or (<X and <Z)")
            .chamfer(2)
            .edges("not %circle")
            .edges("<X or >X or <Y or >Y or <Z or >Z")
            .chamfer(CHAMFER)
        )

        offset_z = REF_HEIGHT - ARM_HOLDER_HEIGHT - 5
        offset_x = REF_OUTER_RADIUS

        self.offset_x = offset_x
        self.offset_z = offset_z
        self.solid = solid


ref = Reflector()
cap = Cap()
arm_holder = ArmHolder()

show_object(ref.solid)
show_object(cap.solid.translate((0, 0, 10)))
show_object(arm_holder.solid.translate((arm_holder.offset_x, 0, arm_holder.offset_z)))

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

ref_export = ref.solid.rotate((0, 0, 0), (0, 1, 0), 45)
ref_export.export(f"{STEP_OUT_DIR}/reflector.step")

cap.solid.export(f"{STEP_OUT_DIR}/cap.step")

holder_export = arm_holder.solid.rotate((0, 0, 0), (0, 1, 0), -90)
holder_export.export(f"{STEP_OUT_DIR}/arm_holder.step")
