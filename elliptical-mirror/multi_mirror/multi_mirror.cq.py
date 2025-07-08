import cadquery as cq
import math


def sub2d(a, b):
    return ((a[0] - b[0]), (a[1] - b[1]))


def hypot2d(a):
    return math.sqrt(a[0] * a[0] + a[1] * a[1])


def add3d(a, b):
    return ((a[0] + b[0]), (a[1] + b[1]), (a[2] + b[2]))


def sub3d(a, b):
    return ((a[0] - b[0]), (a[1] - b[1]), (a[2] - b[2]))


def mul3d(a, s):
    return ((a[0] * s), (a[1] * s), (a[2] * s))


def hypot3d(a):
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def norm3d(a):
    h = hypot3d(a)
    return (a[0] / h, a[1] / h, a[2] / h)


def rotate2d(vec, rad):
    sin = math.sin(rad)
    cos = math.cos(rad)
    x, y = vec
    return (
        x * cos - y * sin,
        x * sin + y * cos,
    )


def cross3d(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def dot3d(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def rotate3d(point, axis_orig, axis_vec, rad):
    va = norm3d(axis_vec)
    vp = sub3d(point, axis_orig)
    cos = math.cos(rad)
    sin = math.sin(rad)
    dot = dot3d(va, vp)
    cross = cross3d(va, vp)
    r_vec = (
        (vp[0] * cos) + (cross[0] * sin) + (va[0] * dot * (1 - cos)),
        (vp[1] * cos) + (cross[1] * sin) + (va[1] * dot * (1 - cos)),
        (vp[2] * cos) + (cross[2] * sin) + (va[2] * dot * (1 - cos)),
    )
    return add3d(axis_orig, r_vec)


def ellipsoid_line_intersection(F1, F2, R_LONG, P1, P2):
    TOL = 1e-8
    T_MIN = -100
    T_MAX = 100
    MAX_ITER = 100

    P_DIST = sub3d(P2, P1)

    def line_func(t):
        return add3d(P1, mul3d(P_DIST, t))

    def dist_func(t):
        q = line_func(t)
        return hypot3d(sub3d(q, F1)) + hypot3d(sub3d(q, F2)) - 2 * R_LONG

    # 線形補間で関数の符号が変わる区間を見つける
    N = 1000
    t_prev = T_MIN
    f_prev = dist_func(t_prev)
    t_cands = []
    for i in range(1, N + 1):
        t_next = T_MIN + (T_MAX - T_MIN) * i / N
        f_next = dist_func(t_next)
        if f_prev * f_next < 0:
            t_cands.append((t_prev, t_next))
        t_prev, f_prev = t_next, f_next

    # この区間に根があるので2分法で解く
    results = []
    for t_prev, t_next in t_cands:
        f_prev = dist_func(t_prev)
        f_next = dist_func(t_next)
        t_best = None
        for i in range(MAX_ITER):
            t_mid = (t_prev + t_next) / 2
            f_mid = dist_func(t_mid)
            t_best = t_mid
            if abs(f_mid) < TOL:
                break
            if f_prev * f_mid < 0:
                t_next, f_next = t_mid, f_mid
            else:
                t_prev, f_prev = t_mid, f_mid
        results.append(line_func(t_best))

    return results


def find_ellipsoid_extremal_point(F1, F2, R_LONG, V):
    TOL = 1e-8
    MAX_ITER = 100

    # V を正規化
    V_norm = hypot3d(V)
    if V_norm == 0:
        raise ValueError("V must not be zero vector!")
    v = [x / V_norm for x in V]

    # 2点F1, F2の中心C
    C = [(F1[i] + F2[i]) * 0.5 for i in range(3)]
    # F1, F2間の距離
    d = hypot3d(sub3d(F2, F1))
    if d > 2 * R_LONG:
        raise ValueError("2a must be greater than |F1-F2|")

    # tについて1次元探索: P = C + t*v
    # t_min, t_maxを推定
    t_min = -2 * R_LONG
    t_max = 2 * R_LONG

    def constraint(t):
        # P = C + t*v
        P = add3d(C, mul3d(v, t))
        return hypot3d(sub3d(P, F1)) + hypot3d(sub3d(P, F2)) - 2 * R_LONG

    # 満たす t を2分法で探索
    # V方向に大きいほう（t_max側）で constraint(t)=0 となる点を探す
    left, right = t_min, t_max
    # まずright側に解があるようにt_maxを広げる
    while constraint(right) < 0:
        right *= 2
        if abs(right) > 1e6:
            raise RuntimeError("Failed to bracket solution on + side")
    # left側に解があるようにt_minを縮める
    while constraint(left) < 0:
        left *= 2
        if abs(left) > 1e6:
            raise RuntimeError("Failed to bracket solution on - side")

    # 2分法
    for i in range(MAX_ITER):
        mid = (left + right) * 0.5
        c = constraint(mid)
        if abs(c) < TOL:
            break
        if c > 0:
            right = mid
        else:
            left = mid
    t_star = (left + right) * 0.5
    P_star = add3d(C, mul3d(v, t_star))
    return P_star


def intersection_plane_line(plane_p, plane_v1, plane_v2, line_p, line_v):
    EPS = 1e-10
    normal = cross3d(plane_v1, plane_v2)
    denom = dot3d(normal, line_v)
    if abs(denom) < EPS:
        return None
    t = dot3d(normal, sub3d(plane_p, line_p)) / denom
    intersection = add3d(line_p, mul3d(line_v, t))
    return intersection


# ----------------------------------------------------------------------

SCALING = 1.5

ENABLE_CHAMFER = True

# 反射板の数
NUM_REFLECTORS = 6

# 反射板の寸法 [mm]
REFLECTOR_INNER_RADIUS = 10
REFLECTOR_OUTER_RADIUS = 60 * SCALING
REFLECTOR_TOP_CUT_SIZE = 10 * SCALING
REFLECTOR_THICKNESS = 3

# 反射板のサポート有無
REFLECTOR_ENABLE_SUPPORT = True

# 焦点間の距離 [mm]
FOCUS_DISTANCE = 40 * SCALING

# 楕円の傾き [度]
ELLIPSE_ANGLE = 45

# 楕円の半径 [mm]
ELLIPSE_SHORT_RADIUS = 40 * SCALING
ELLIPSE_LONG_RADIUS = math.sqrt((FOCUS_DISTANCE / 2) ** 2 + ELLIPSE_SHORT_RADIUS**2)

REFLECTOR_FLOOR_Z = -10 * SCALING

# 焦点の座標
FOCUS_NEAR = rotate3d(
    (FOCUS_DISTANCE, 0, 0), (0, 0, 0), (0, -1, 0), math.radians(ELLIPSE_ANGLE)
)
FOCUS_FAR = rotate3d((0, 0, 0), (0, 0, 0), (0, -1, 0), math.radians(ELLIPSE_ANGLE))

# フレームの天井の高さ
ellipse_top = find_ellipsoid_extremal_point(
    FOCUS_NEAR, FOCUS_FAR, ELLIPSE_LONG_RADIUS, (0, 0, 1)
)
FRAME_CEIL_Z = math.ceil(ellipse_top[2] + 10)
FRAME_THICKNESS = 3

# 嵌合用の隙間 [mm]
GENERIC_GAP = 0.5

# 面取りのサイズ [mm]
CHAMFER = 0.5

# 天井部分の固定用穴の間隔 [mm]
CEIL_INNER_MOUNT_X = 15
CEIL_OUTER_MOUNT_X = 30 * SCALING
CEIL_OUTER_MOUNT_Y = 20
CEIL_EXTRA_MOUNT_X = REFLECTOR_OUTER_RADIUS - REFLECTOR_TOP_CUT_SIZE - 5

# 固定用キャップの直径 [mm]
CAP_DIAMETER = (CEIL_OUTER_MOUNT_X + 15) * 2

# 固定用キャップの穴の間隔 [mm]
CAP_MOUNT_DISTANCE_X = 50
CAP_MOUNT_DISTANCE_Y = 30

# LED アームの取り付け用穴の間隔 [mm]
ARM_MOUNT_HOLE_DISTANCE = 10

ARM_BASE_HEIGHT = 20
ARM_MOUNT_Z_OFFSET = 10 * SCALING

ARM_JOINT_DISTANCE = 25

ARM_LINK_WIDTH = 12
ARM_LINK_LENGTH = ARM_JOINT_DISTANCE + ARM_LINK_WIDTH
ARM_LINK_THICKNESS = 10

# ボールジョイントの寸法
BALLJOINT_NECK_THICKNESS = 5
BALLJOINT_NECK_WIDTH = 20
BALLJOINT_NECK_HEIGHT = 3
BALLJOINT_BALL_DIAMETER = BALLJOINT_NECK_WIDTH * 1.25
BALLJOINT_CUT_DIAMETER = 15
BALLJOINT_HOLE_DIAMETER = 5

M3_TAP_HOLE_DIAMETER = 2.5
M3_HOLE_DIAMETER = 3 + GENERIC_GAP
M3_HEAD_DIAMETER = 5.5
M3_HEAD_THICKNESS = 2.5
M3_NUT_DIAMETER = 5.5 * 2 / math.sqrt(3) + 0.2
M3_NUT_THICKNESS = 2.5
M3_BUTTERFLY_NUT_THICKNESS = 3.5
M3_WASHER_THICKNESS = 0.5

STEP_OUT_DIR = "./step"


class MirrorSegment:
    def __init__(self):
        inner_x = REFLECTOR_INNER_RADIUS * math.cos(math.pi / NUM_REFLECTORS)
        inner_y = REFLECTOR_INNER_RADIUS * math.sin(math.pi / NUM_REFLECTORS)
        solid = (
            cq.Workplane("XY")
            .ellipseArc(
                REFLECTOR_OUTER_RADIUS,
                REFLECTOR_OUTER_RADIUS,
                -180 / NUM_REFLECTORS,
                180 / NUM_REFLECTORS,
                startAtCurrent=False,
            )
            .lineTo(inner_x, inner_y)
            .lineTo(inner_x, -inner_y)
            .close()
            .extrude(FRAME_CEIL_Z)
        )

        # 楕円面を彫り込む
        cutter = (
            cq.Workplane("XY")
            .box(400, 400, 200, centered=(True, True, False))
            .translate((0, 0, -200 + REFLECTOR_FLOOR_Z))
        )
        cutter = cutter.union(
            cq.Workplane("XY")
            .ellipseArc(
                ELLIPSE_LONG_RADIUS,
                ELLIPSE_SHORT_RADIUS,
                0,
                180,
                sense=-1,
                startAtCurrent=False,
            )
            .close()
            .revolve(360, (0, 0, 0), (1, 0, 0))
            .translate((FOCUS_DISTANCE / 2, 0, 0))
        )
        cutter = cutter.rotate((0, 0, 0), (0, -1, 0), ELLIPSE_ANGLE)
        solid = solid.cut(cutter)

        # LED 取り付け用の溝を掘る
        y_pos = REFLECTOR_OUTER_RADIUS - 5
        z = 10
        verts = [
            (y_pos, -z),
            (y_pos, z),
            (y_pos + 10, z + 10),
            (y_pos + 10, -z - 10),
        ]
        cutter = cq.Workplane("XY").polyline(verts).close().extrude(FRAME_CEIL_Z + 10)
        solid = solid.cut(cutter)

        # 天井を斜めにする
        verts = [
            (REFLECTOR_OUTER_RADIUS - REFLECTOR_TOP_CUT_SIZE, FRAME_CEIL_Z + 10),
            (REFLECTOR_OUTER_RADIUS, FRAME_CEIL_Z + 10),
            (REFLECTOR_OUTER_RADIUS, FRAME_CEIL_Z - REFLECTOR_TOP_CUT_SIZE),
            (REFLECTOR_OUTER_RADIUS - REFLECTOR_TOP_CUT_SIZE, FRAME_CEIL_Z),
        ]
        cutter = (
            cq.Workplane("XZ")
            .polyline(verts)
            .close()
            .extrude(REFLECTOR_OUTER_RADIUS, both=True)
        )
        solid = solid.cut(cutter)

        # 穴開け
        solid = (
            solid
            # キャップ用の穴
            .faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (CEIL_INNER_MOUNT_X, 0),
                    (CEIL_OUTER_MOUNT_X, -CEIL_OUTER_MOUNT_Y),
                    (CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y),
                    (CEIL_EXTRA_MOUNT_X, 0),
                ]
            )
            .circle(M3_TAP_HOLE_DIAMETER / 2)
            .cutBlind(-5)
            # LED アーム用の穴
            .faces(">>X[-2]")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (0, FRAME_CEIL_Z - ARM_MOUNT_Z_OFFSET),
                    (
                        0,
                        FRAME_CEIL_Z - ARM_MOUNT_Z_OFFSET - ARM_MOUNT_HOLE_DISTANCE,
                    ),
                ]
            )
            .circle(M3_TAP_HOLE_DIAMETER / 2)
            .cutBlind(-5)
        )

        # 面取り
        if False and ENABLE_CHAMFER:
            solid = (
                solid.faces("<X or <<Y[-2] or >>Y[-2] or <<Z[-2]")
                # .faces("<X or <<Y[-2] or >>Y[-2] or <<Z[-3]")
                .chamfer(CHAMFER)
                # .edges("%circle")
                # .chamfer(CHAMFER)
            )

        intersector = (
            cq.Workplane("XY")
            .ellipseArc(
                ELLIPSE_LONG_RADIUS + REFLECTOR_THICKNESS,
                ELLIPSE_SHORT_RADIUS + REFLECTOR_THICKNESS,
                0,
                180,
                sense=-1,
                startAtCurrent=False,
            )
            .close()
            .revolve(360, (0, 0, 0), (1, 0, 0))
            .translate((FOCUS_DISTANCE / 2, 0, 0))
            .rotate((0, 0, 0), (0, -1, 0), ELLIPSE_ANGLE)
        )

        # 反射板を切り出す
        reflector = solid.intersect(intersector)

        reflector_support = None

        # サポートを反射面に沿って切り取るための楕円のパラメータ
        f1 = (-FOCUS_DISTANCE / 2, 0, 0)
        f2 = (FOCUS_DISTANCE / 2, 0, 0)
        r_short = ELLIPSE_SHORT_RADIUS + REFLECTOR_THICKNESS - 0.1
        r_long = ELLIPSE_LONG_RADIUS + REFLECTOR_THICKNESS - 0.1

        # サポートの取り付け角度
        angle_deg_list = [0, 20, 40]

        # サポートの高さ
        height_list = [r_short * 3 / 4, r_short * 2 / 3, r_short / 2]

        for angle_deg, height in zip(angle_deg_list, height_list):
            angle_rad = math.radians(angle_deg)
            cos = math.cos(angle_rad)
            sin = math.sin(angle_rad)
            top_z = REFLECTOR_FLOOR_Z + height

            # サポートの頂点位置の算出
            p1 = add3d(f2, (0, 0, top_z))
            p2 = add3d(p1, (r_short * cos, r_short * sin, 0))
            intersections = ellipsoid_line_intersection(f1, f2, r_long, p1, p2)
            top = max(intersections, key=lambda p: p[0])

            # スケッチのため XZ 平面に回転
            top = rotate3d(top, p1, (0, 0, 1), -angle_rad)

            # サポートのソリッド生成
            foot_x = top[0] + height * 2 / 3
            foot_diameter = 15
            verts = [
                (p1[0], REFLECTOR_FLOOR_Z),
                (top[0], top_z),
                (foot_x, REFLECTOR_FLOOR_Z),
            ]
            support = (
                cq.Workplane("XZ").polyline(verts).close().extrude(0.25, both=True)
            )
            support = support.union(
                cq.Workplane("XY")
                .pushPoints([(foot_x - foot_diameter / 2, 0)])
                .circle(foot_diameter / 2)
                .extrude(1)
                .translate((0, 0, REFLECTOR_FLOOR_Z))
            )

            if False:
                # ミシン目を付ける
                hole_diameter = 2
                num_holes = int(math.floor(height / (hole_diameter * 2))) - 3
                for i in range(num_holes):
                    z = REFLECTOR_FLOOR_Z + (i + 1) * (hole_diameter * 2)
                    p1 = add3d(f2, (0, 0, z))
                    p2 = add3d(p1, (r_short * cos, r_short * sin, 0))
                    intersections = ellipsoid_line_intersection(f1, f2, r_long, p1, p2)
                    pos = max(intersections, key=lambda p: p[0])
                    pos = rotate3d(pos, p1, (0, 0, 1), -angle_rad)
                    support = support.cut(
                        cq.Workplane("XY")
                        .box(hole_diameter, 10, hole_diameter)
                        .translate(pos)
                    )

            # 取り付け位置に移動
            support = support.rotate(f2, add3d(f2, (0, 0, 1)), angle_deg).translate(
                (FOCUS_DISTANCE / 2, 0, 0)
            )

            # 楕円面に沿って切り取る
            support = support.cut(
                cq.Workplane("XY")
                .ellipseArc(
                    r_long,
                    r_short,
                    0,
                    180,
                    sense=-1,
                    startAtCurrent=False,
                )
                .close()
                .revolve(360, (0, 0, 0), (1, 0, 0))
                .translate((FOCUS_DISTANCE / 2, 0, 0))
            )

            # 反射板の位置に合わせる
            support = support.rotate((0, 0, 0), (0, -1, 0), ELLIPSE_ANGLE)

            # 左右対称にする
            if angle_deg != 0:
                support = support.union(support.mirror("XZ"))

            if reflector_support is None:
                reflector_support = support
            else:
                reflector_support = reflector_support.union(support)

        # フレームを切り出す
        frame = solid.cut(intersector.translate((GENERIC_GAP, 0, GENERIC_GAP)))

        outer_x = REFLECTOR_OUTER_RADIUS * math.cos(math.pi / NUM_REFLECTORS)
        outer_y = REFLECTOR_OUTER_RADIUS * math.sin(math.pi / NUM_REFLECTORS)
        verts = [
            (inner_x, -inner_y),
            (inner_x, inner_y),
            (outer_x, outer_y),
            (outer_x, -outer_y),
        ]

        cutter = cq.Workplane("XY").polyline(verts).close().extrude(FRAME_CEIL_Z + 1)
        cutter = cutter.cut(
            cq.Workplane("XY")
            .box(400, 400, 400, centered=(False, True, False))
            .translate((CEIL_INNER_MOUNT_X - 400 + 5, 0, 0))
        )
        cutter = cutter.cut(
            cq.Workplane("XY")
            .box(400, 400, 400, centered=(False, True, False))
            .translate((CEIL_EXTRA_MOUNT_X, 0, 0))
        )

        # 固定用穴の周辺
        hole_mount = (
            cq.Workplane("XY")
            .box(6, 100, FRAME_CEIL_Z, centered=(True, False, False))
            .translate((CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y, 0))
            .faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y)])
            .circle(6 / 2)
            .extrude(-100)
        )
        cutter = cutter.cut(hole_mount)
        cutter = cutter.cut(hole_mount.mirror("XZ"))

        # トラス
        column_height = 5
        truss = (
            cq.Workplane("XY")
            .box(
                FRAME_THICKNESS,
                CEIL_OUTER_MOUNT_Y * 2,
                column_height,
                centered=(True, True, False),
            )
            .translate((CEIL_OUTER_MOUNT_X, 0, FRAME_CEIL_Z - column_height))
        )
        p1 = (CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y, FRAME_CEIL_Z - column_height)
        p2 = (CEIL_EXTRA_MOUNT_X, 0, FRAME_CEIL_Z - column_height)
        column_angle_deg = math.degrees(math.atan2(p2[1] - p1[1], p2[0] - p1[0])) - 90
        column_length = hypot3d(sub3d(p2, p1))
        column = (
            cq.Workplane("XY")
            .box(
                FRAME_THICKNESS,
                column_length,
                column_height,
                centered=(True, False, False),
            )
            .rotate((0, 0, 0), (0, 0, 1), column_angle_deg)
            .translate(p1)
        )
        truss = truss.union(column)
        truss = truss.union(column.mirror("XZ"))
        truss = truss.union(
            cq.Workplane("XY")
            .box(10, 10, column_height, centered=(True, True, False))
            .translate(p2)
        )
        cutter = cutter.cut(truss)

        # 反射面とフレームの間に隙間を空ける
        ellip = (
            cq.Workplane("XY")
            .ellipseArc(
                ELLIPSE_LONG_RADIUS + REFLECTOR_THICKNESS + 5,
                ELLIPSE_SHORT_RADIUS + REFLECTOR_THICKNESS + 5,
                0,
                180,
                sense=-1,
                startAtCurrent=False,
            )
            .close()
            .revolve(360, (0, 0, 0), (1, 0, 0))
            .translate((FOCUS_DISTANCE / 2, 0, 0))
        )
        ellip = ellip.cut(
            cq.Workplane("XY")
            .box(400, 400, 200, centered=(True, True, False))
            .translate((0, 0, -200 + REFLECTOR_FLOOR_Z + 5))
        )
        ellip = ellip.rotate((0, 0, 0), (0, -1, 0), ELLIPSE_ANGLE)
        ellip = ellip.cut(
            cq.Workplane("XY")
            .box(400, 400, 400, centered=(False, True, True))
            .translate((CEIL_EXTRA_MOUNT_X - 400 - 1, 0, 0))
        )
        cutter = cutter.union(ellip)

        # 側面のフレーム部分を残す
        cutter = cutter.cut(
            cq.Workplane("XY")
            .box(400, 400, 400, centered=(False, False, False))
            .translate((0, -FRAME_THICKNESS, 0))
            .rotate((0, 0, 0), (0, 0, 1), 180 / NUM_REFLECTORS)
        )
        cutter = cutter.cut(
            cq.Workplane("XY")
            .box(400, 400, 400, centered=(False, False, False))
            .translate((0, -400 + FRAME_THICKNESS, 0))
            .rotate((0, 0, 0), (0, 0, 1), -180 / NUM_REFLECTORS)
        )

        # フレームの穴部分を削る
        cutter = cutter.edges("|Z").edges(">Y or <Y").chamfer(2)
        frame = frame.cut(cutter)

        # 側面の穴
        x = CEIL_EXTRA_MOUNT_X
        z = FRAME_CEIL_Z - 5
        verts = [
            (x, z),
            (x - 20, z),
            (x - 20, z - 5),
            (x, z - 25),
        ]
        cutter = (
            cq.Workplane("XZ")
            .polyline(verts)
            .close()
            .extrude(REFLECTOR_OUTER_RADIUS, both=True)
            .edges("|Y")
            .edges(">X or >Z")
            .chamfer(2)
        )
        frame = frame.cut(cutter)

        # 外周の穴
        verts = [
            (0, 0),
            (25, -8),
            (25, -100),
            (0, -100),
        ]
        hole = (
            cq.Workplane("YZ")
            .polyline(verts)
            .close()
            .extrude(100)
            .translate(
                (
                    REFLECTOR_OUTER_RADIUS - REFLECTOR_TOP_CUT_SIZE - 10,
                    FRAME_THICKNESS,
                    FRAME_CEIL_Z - 10,
                )
            )
            .edges("|X")
            .chamfer(2)
            .rotate((0, 0, 0), (0, 0, 1), -180 / NUM_REFLECTORS)
        )
        frame = frame.cut(hole)
        frame = frame.cut(hole.mirror("XZ"))

        self.reflector_solid = reflector
        self.reflector_support = reflector_support
        self.frame_solid = frame
        self.bounding_solid = solid


class HeadParams:
    NORMAL = 0
    DISH = 1
    CAP = 2
    NUT = 3

    def __init__(self, type=NORMAL, diameter=0, thickness=0):
        self.type = type
        self.diameter = diameter
        self.thickness = thickness

        if type == HeadParams.DISH:
            verts = [
                (0, 0),
                (diameter / 2, 0),
                (0, -diameter / 2),
            ]
            self.solid = (
                cq.Workplane("XZ")
                .polyline(verts)
                .close()
                .revolve(360, (0, 0, 0), (0, 1, 0))
            )
        elif type == HeadParams.CAP:
            verts = [
                (0, 0),
                (diameter / 2, 0),
                (diameter / 2, -thickness),
                (0, -thickness),
            ]
            self.solid = (
                cq.Workplane("XZ")
                .polyline(verts)
                .close()
                .revolve(360, (0, 0, 0), (0, 1, 0))
            )
        elif type == HeadParams.NUT:
            self.solid = (
                cq.Workplane("XY")
                .polygon(6, diameter * 2 / math.sqrt(3))
                .extrude(-thickness)
            )
        else:
            self.solid = None


def create_hole_cutter(
    diameter,
    length,
    near: HeadParams = None,
    far: HeadParams = None,
):
    solid = (
        cq.Workplane("XY")
        .cylinder(length, diameter / 2, centered=(True, True, False))
        .translate((0, 0, -length))
    )

    if near is not None and near.type != HeadParams.NORMAL:
        solid = solid.union(near.solid)

    if far is not None and far.type != HeadParams.NORMAL:
        solid = solid.union(
            far.solid.rotate((0, 0, 0), (1, 0, 0), 180).translate((0, 0, -length))
        )

    return solid


M3_Dish = HeadParams(HeadParams.DISH, diameter=7)
M3_Nut = HeadParams(HeadParams.NUT, diameter=5.5, thickness=2.5)
M3_Chamfer = HeadParams(HeadParams.DISH, diameter=M3_HOLE_DIAMETER + CHAMFER * 2)


class MirrorFastener:
    def __init__(self):
        center = (-1, 0, FRAME_CEIL_Z)

        beam_thickness = 5
        beam_width_narrow = 10
        beam_width_wide = 20
        beam_length = REFLECTOR_OUTER_RADIUS - 10

        hole_points = [
            (-1, -beam_width_wide / 2),
            (-1, beam_width_wide / 2),
            (beam_length - 2, beam_width_narrow / 2),
            (beam_length, 0),
            (beam_length - 2, -beam_width_narrow / 2),
        ]
        beam_template = (
            cq.Workplane("XY")
            .polyline(hole_points)
            .close()
            .extrude(beam_thickness)
            .translate((0, 0, FRAME_CEIL_Z))
        )

        hole_points = [
            (beam_length - 20, beam_thickness + 0.1),
            (beam_length + 0.1, beam_thickness + 0.1),
            (beam_length + 0.1, 2),
        ]
        beam_template = beam_template.cut(
            cq.Workplane("XZ")
            .polyline(hole_points)
            .close()
            .extrude(beam_width_narrow, both=True)
            .translate((0, 0, FRAME_CEIL_Z))
        )

        beam_template = beam_template.rotate(
            center, add3d(center, (0, 0, 1)), 180 / NUM_REFLECTORS
        )

        solid = beam_template
        for i in range(1, NUM_REFLECTORS):
            bam_rotated = beam_template.rotate(
                center, add3d(center, (0, 0, 1)), i * 360 / NUM_REFLECTORS
            )
            if solid is None:
                solid = bam_rotated
            else:
                solid = solid.union(bam_rotated)

        # 中心の板と穴
        d = REFLECTOR_INNER_RADIUS * 2 * math.sqrt(3) / 2
        solid = (
            solid.faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .polygon(NUM_REFLECTORS, d + 15 * 2, circumscribed=True)
            .extrude(-beam_thickness)
            .faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .polygon(NUM_REFLECTORS, d, circumscribed=True)
            .cutBlind(-100)
        )

        # ねじ穴のテンプレート
        hole_cutter = create_hole_cutter(
            M3_HOLE_DIAMETER,
            length=beam_thickness,
            near=M3_Dish,
            far=M3_Chamfer,
        )

        # 反射板用のネジ穴開け
        hole_points = [
            (CEIL_INNER_MOUNT_X, 0, FRAME_CEIL_Z + beam_thickness),
            (CEIL_OUTER_MOUNT_X, -CEIL_OUTER_MOUNT_Y, FRAME_CEIL_Z + beam_thickness),
            (CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y, FRAME_CEIL_Z + beam_thickness),
        ]
        for i in range(NUM_REFLECTORS):
            for p in hole_points:
                p_rotated = rotate3d(
                    p, center, (0, 0, 1), math.radians(i * 360 / NUM_REFLECTORS)
                )
                solid = solid.union(
                    cq.Workplane("XY")
                    .cylinder(beam_thickness, 5, centered=(True, True, False))
                    .translate(p_rotated)
                    .translate((0, 0, -beam_thickness))
                )
                solid = solid.cut(hole_cutter.translate(p_rotated))

        # ねじ穴のテンプレート
        hole_cutter = create_hole_cutter(
            M3_HOLE_DIAMETER,
            length=beam_thickness,
            near=M3_Dish,
            far=M3_Chamfer,
        ).rotate((0, 0, 0), (1, 0, 0), 180)

        # マウント用の穴開け
        hole_points = [
            add3d(center, (CAP_MOUNT_DISTANCE_X / 2, -CAP_MOUNT_DISTANCE_Y / 2, 0)),
            add3d(center, (CAP_MOUNT_DISTANCE_X / 2, CAP_MOUNT_DISTANCE_Y / 2, 0)),
            add3d(center, (-CAP_MOUNT_DISTANCE_X / 2, CAP_MOUNT_DISTANCE_Y / 2, 0)),
            add3d(center, (-CAP_MOUNT_DISTANCE_X / 2, -CAP_MOUNT_DISTANCE_Y / 2, 0)),
        ]
        for p in hole_points:
            solid = solid.cut(hole_cutter.translate(p))

        self.solid = solid


def create_ball_joint():
    solid = cq.Workplane("XY").sphere(BALLJOINT_BALL_DIAMETER / 2)
    solid = solid.intersect(
        cq.Workplane("XY")
        .cylinder(BALLJOINT_BALL_DIAMETER * 2, BALLJOINT_CUT_DIAMETER / 2)
        .rotate((0, 0, 0), (1, 0, 0), 90)
    )
    solid = solid.cut(
        cq.Workplane("XY")
        .cylinder(BALLJOINT_BALL_DIAMETER * 2, BALLJOINT_HOLE_DIAMETER / 2)
        .rotate((0, 0, 0), (1, 0, 0), 90)
    )
    solid = solid.translate((0, 0, BALLJOINT_NECK_HEIGHT + BALLJOINT_CUT_DIAMETER / 2))
    solid = solid.union(
        cq.Workplane("XY").box(
            BALLJOINT_NECK_THICKNESS,
            BALLJOINT_NECK_WIDTH,
            BALLJOINT_NECK_HEIGHT + 5,
            centered=(True, True, False),
        )
    )
    return solid


class LedArm:
    def __init__(self):
        # ベース部分
        base_thickness = 5

        joint1_width = 20
        joint1_height = ARM_MOUNT_HOLE_DISTANCE + 5 * 2 + 10

        joint1_solid = (
            cq.Workplane("XY")
            .box(
                base_thickness,
                joint1_width,
                joint1_height,
                centered=(False, True, False),
            )
            .faces(">X")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (0, 5),
                    (0, 5 + ARM_MOUNT_HOLE_DISTANCE),
                ]
            )
            .circle(M3_HOLE_DIAMETER / 2)
            .cutBlind(-10)
            # 取り付け穴の皿ネジ用の座繰り
            .edges("%circle and >X")
            .chamfer(2)
        )

        # グリップを生成してくっつける
        joint1_solid = joint1_solid.union(
            create_ball_joint()
            .rotate((0, 0, 0), (0, 1, 0), 90)
            .translate(
                (base_thickness, 0, joint1_height - BALLJOINT_NECK_THICKNESS / 2)
            )
        )

        # 全体の面取り
        if ENABLE_CHAMFER:
            joint1_solid = joint1_solid.faces("-X or +Y or -Y or -Z or +Z").chamfer(
                CHAMFER
            )

        ARM2_WIDTH = 10
        HEAT_SINK_HOLE_DISTANCE = 30

        hole_points = [
            (0, HEAT_SINK_HOLE_DISTANCE / 2),
            (0, -HEAT_SINK_HOLE_DISTANCE / 2),
        ]

        joint2_solid = (
            cq.Workplane("XY")
            .box(
                ARM2_WIDTH,
                HEAT_SINK_HOLE_DISTANCE,
                base_thickness,
                centered=(True, True, False),
            )
            # 両端の丸み
            .pushPoints(hole_points)
            .circle(ARM2_WIDTH / 2)
            .extrude(base_thickness)
            # ヒートシンク取り付け穴
            .faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints(hole_points)
            .circle(M3_HOLE_DIAMETER / 2)
            .cutBlind(-100)
            # ナット用の凹み
            .pushPoints(hole_points)
            .polygon(6, M3_NUT_DIAMETER)
            .cutBlind(-M3_NUT_THICKNESS)
        )
        joint2_solid = joint2_solid.union(
            create_ball_joint().translate((0, 0, base_thickness))
        )
        joint2_solid = joint2_solid.translate(
            (
                base_thickness + BALLJOINT_NECK_HEIGHT + BALLJOINT_CUT_DIAMETER / 2,
                0,
                joint1_height
                - ARM_JOINT_DISTANCE
                - BALLJOINT_CUT_DIAMETER
                - BALLJOINT_NECK_HEIGHT,
            )
        )

        # 全体の面取り
        if ENABLE_CHAMFER:
            joint2_solid = joint2_solid.edges("%line and (<X or >X)").chamfer(CHAMFER)

        link_solid = (
            cq.Workplane("XY")
            .box(
                ARM_LINK_WIDTH,
                ARM_LINK_THICKNESS,
                ARM_LINK_LENGTH - ARM_LINK_WIDTH,
                centered=(True, False, True),
            )
            # 両端の丸み
            .faces(">Y")
            .workplane(origin=(0, 0, 0))
            .pushPoints(
                [
                    (0, ARM_JOINT_DISTANCE / 2),
                    (0, -ARM_JOINT_DISTANCE / 2),
                ]
            )
            .circle(ARM_LINK_WIDTH / 2)
            .extrude(-ARM_LINK_THICKNESS)
            # ねじ穴部分の構造強化
            .faces(">Y")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .circle(ARM_LINK_WIDTH / 2 + 2)
            .extrude(-ARM_LINK_THICKNESS)
            # 中央の穴
            .faces(">Y")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .circle(M3_HOLE_DIAMETER / 2)
            .cutBlind(-100)
        )
        # ボールジョイントの形に削る
        cutter = (
            cq.Workplane("XY")
            .sphere(BALLJOINT_BALL_DIAMETER / 2)
            .translate((0, -BALLJOINT_NECK_WIDTH / 2, 0))
        )
        cutter = cutter.cut(
            cq.Workplane("XY").cylinder(100, 1).rotate((0, 0, 0), (1, 0, 0), 90)
        )
        cutter = cutter.translate((0, 0, ARM_JOINT_DISTANCE / 2))
        link_solid = link_solid.cut(cutter).cut(cutter.mirror("XY"))

        left_link_solid = (
            link_solid
            # ねじ頭が埋まる穴
            .faces(">Y")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .circle((M3_HEAD_DIAMETER + GENERIC_GAP) / 2)
            .cutBlind(-M3_BUTTERFLY_NUT_THICKNESS - M3_WASHER_THICKNESS - GENERIC_GAP)
            # ナットがはまる壁
            .faces("<Y")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .polygon(6, M3_NUT_DIAMETER + 3)
            .extrude(M3_NUT_THICKNESS)
            .faces("<Y")
            .workplane(origin=(0, 0, 0))
            .pushPoints([(0, 0)])
            .polygon(6, M3_NUT_DIAMETER)
            .cutBlind(-M3_NUT_THICKNESS)
        )
        right_link_solid = link_solid

        # 全体の面取り
        if ENABLE_CHAMFER:
            left_link_solid = (
                left_link_solid.edges("%circle and >>Y")
                .edges("not (>>Z or <<Z or >>X or <<X)")
                .chamfer(CHAMFER)
                .edges(">>Y and (>>Z or <<Z or >>X or <<X)")
                .chamfer(2)
            )
            right_link_solid = (
                right_link_solid.edges("%circle and >>Y")
                .edges("not (>>Z or <<Z or >>X or <<X)")
                .chamfer(CHAMFER)
                .edges(">>Y and (>>Z or <<Z or >>X or <<X)")
                .chamfer(2)
            )

        offset = (
            base_thickness + BALLJOINT_NECK_HEIGHT + BALLJOINT_CUT_DIAMETER / 2,
            BALLJOINT_NECK_WIDTH / 2,
            joint1_height - (ARM_JOINT_DISTANCE + BALLJOINT_CUT_DIAMETER) / 2 + 5,
        )
        left_link_solid = left_link_solid.translate(offset)
        right_link_solid = right_link_solid.translate(offset).mirror("XZ")

        # 然るべき位置に移動
        offset = (
            REFLECTOR_OUTER_RADIUS - 5,
            0,
            FRAME_CEIL_Z - ARM_BASE_HEIGHT - ARM_MOUNT_Z_OFFSET + 5,
        )
        joint1_solid = joint1_solid.translate(offset)
        joint2_solid = joint2_solid.translate(offset)
        left_link_solid = left_link_solid.translate(offset)
        right_link_solid = right_link_solid.translate(offset)

        self.joint1_solid = joint1_solid
        self.joint2_solid = joint2_solid
        self.left_link_solid = left_link_solid
        self.right_link_solid = right_link_solid


class FocusIndicator:
    def __init__(self, body: MirrorSegment):
        # 内側の足の座標
        # (inner_edge_x, inner_edge_y) = rotate2d((REFLECTOR_INNER_RADIUS, 0), math.pi / NUM_REFLECTORS)
        (inner_edge_x, inner_edge_y) = (
            REFLECTOR_INNER_RADIUS * math.cos(math.pi / NUM_REFLECTORS),
            0,
        )
        inner_cross = ellipsoid_line_intersection(
            FOCUS_NEAR,
            FOCUS_FAR,
            ELLIPSE_LONG_RADIUS,
            (inner_edge_x, inner_edge_y, 0),
            (inner_edge_x, inner_edge_y, FRAME_CEIL_Z),
        )
        foot_pos1 = max(inner_cross, key=lambda p: p[2])

        # 外側のエッジ
        outer_edge_p = rotate3d(
            (REFLECTOR_OUTER_RADIUS, 0, 0),
            (0, 0, 0),
            (0, 0, 1),
            math.pi / NUM_REFLECTORS,
        )
        outer_edge_v = (0, 0, 1)

        # 反射板の壁の下面の切り取る面
        plane_rad = math.radians(ELLIPSE_ANGLE)
        plane_p1 = rotate3d((0, 0, REFLECTOR_FLOOR_Z), (0, 0, 0), (0, -1, 0), plane_rad)
        plane_p2 = rotate3d((1, 0, REFLECTOR_FLOOR_Z), (0, 0, 0), (0, -1, 0), plane_rad)
        plane_p = plane_p1
        plane_v1 = sub3d(plane_p2, plane_p1)
        plane_v2 = (0, 1, 0)

        # エッジと面の交点
        intersector_p1 = plane_p
        intersector_p2 = intersection_plane_line(
            plane_p,
            plane_v1,
            plane_v2,
            outer_edge_p,
            outer_edge_v,
        )

        # 楕円体との交点
        intersections = ellipsoid_line_intersection(
            FOCUS_NEAR,
            FOCUS_FAR,
            ELLIPSE_LONG_RADIUS,
            intersector_p1,
            intersector_p2,
        )
        foot_pos2 = max(intersections, key=lambda p: p[2])
        foot_pos2 = add3d(
            foot_pos2, mul3d(norm3d(sub3d(foot_pos1, foot_pos2)), 2 * SCALING)
        )

        # スケッチのため YZ 平面に回転
        sketch_angle_rad = math.atan2(
            foot_pos2[2] - foot_pos1[2], foot_pos2[0] - foot_pos1[0]
        )
        sketch_angle_rad += math.pi / 2
        sketch_angle_deg = math.degrees(sketch_angle_rad)
        foot_pos1_rot = rotate3d(foot_pos1, (0, 0, 0), (0, -1, 0), -sketch_angle_rad)
        foot_pos2_rot = rotate3d(foot_pos2, (0, 0, 0), (0, -1, 0), -sketch_angle_rad)
        foot_pos3_rot = (foot_pos2_rot[0], -foot_pos2_rot[1], foot_pos2_rot[2])
        near_rot = rotate3d(FOCUS_NEAR, (0, 0, 0), (0, -1, 0), -sketch_angle_rad)
        far_rot = rotate3d(FOCUS_FAR, (0, 0, 0), (0, -1, 0), -sketch_angle_rad)

        wall_height = 10
        wall_thickness = 2

        # 足を形成
        foot_radius = 5 * SCALING
        foot1_pos_xy = (foot_pos1_rot[1], foot_pos1_rot[2])
        foot2_pos_xy = (foot_pos2_rot[1], foot_pos2_rot[2])
        foot3_pos_xy = (-foot_pos2_rot[1], foot_pos2_rot[2])
        foots = (
            cq.Workplane("YZ")
            .pushPoints([foot1_pos_xy, foot2_pos_xy, foot3_pos_xy])
            .circle(foot_radius)
            .extrude(wall_height / 2, both=True)
            .translate((foot_pos2_rot[0], 0, 0))
        )

        wall12_length = hypot3d(sub3d(foot_pos2_rot, foot_pos1_rot))
        wall12_angle = math.atan2(
            foot_pos2_rot[2] - foot_pos1_rot[2], foot_pos2_rot[1] - foot_pos1_rot[1]
        )
        wall12_angle = math.degrees(wall12_angle) + 90
        wall1 = (
            cq.Workplane("ZX")
            .rect(wall12_length, wall_height, centered=(False, True))
            .extrude(wall_thickness / 2, both=True)
            .translate((foot_pos1_rot[0], 0, foot_pos1_rot[2] - wall12_length))
            .rotate(
                (0, foot_pos1_rot[1], foot_pos1_rot[2]),
                (1, foot_pos1_rot[1], foot_pos1_rot[2]),
                wall12_angle,
            )
        )
        wall2 = wall1.mirror("XZ")

        wall3_length = (foot_pos3_rot[1] - foot_pos2_rot[1]) / 2
        wall3_z = foot_pos2_rot[2] + (foot_pos1_rot[2] - foot_pos2_rot[2]) / 2
        wall3 = (
            cq.Workplane("YX")
            .rect(wall3_length, wall_height, centered=(True, True))
            .extrude(wall_thickness / 2, both=True)
            .translate((foot_pos2_rot[0], 0, wall3_z))
        )

        arrow_thickness = 1
        arrow_size = 2 * SCALING
        arrow_space_size = arrow_size + 5 * SCALING
        arrow_body_size = arrow_space_size + 3 * SCALING
        verts = [
            (arrow_body_size, -arrow_size),
            (arrow_space_size, -arrow_size),
            (-arrow_size, arrow_space_size),
            (-arrow_size, arrow_body_size),
            (arrow_body_size, arrow_body_size),
        ]
        arrow_body = (
            cq.Workplane("XZ")
            .polyline(verts)
            .close()
            .extrude(arrow_thickness / 2, both=True)
        )
        verts = [
            (arrow_space_size - arrow_size * 2, -arrow_size),
            (arrow_space_size, arrow_size),
            (arrow_space_size, arrow_space_size),
            (arrow_size, arrow_space_size),
            (-arrow_size, arrow_space_size - arrow_size * 2),
        ]
        arrow_cutter = (
            cq.Workplane("XZ")
            .polyline(verts)
            .close()
            .extrude(wall_thickness, both=True)
        )

        floor_x = foot_pos1_rot[0] + wall_height / 2
        foot1_z = foot_pos1_rot[2]
        (near_x, _, near_z) = near_rot
        (far_x, _, far_z) = far_rot
        verts = [
            (floor_x, foot1_z),
            (floor_x, near_z),
            (near_x + arrow_body_size, near_z - arrow_size),
            (near_x + arrow_space_size, near_z + arrow_space_size),
            (near_x - arrow_size, near_z + arrow_body_size),
            (far_x + arrow_body_size, far_z - arrow_size),
            (far_x + arrow_space_size, far_z + arrow_space_size),
            (far_x - arrow_size, far_z + arrow_body_size),
        ]
        indicator = (
            cq.Workplane("XZ")
            .polyline(verts)
            .close()
            .extrude(wall_thickness / 2, both=True)
        )

        far_body = arrow_body.translate(far_rot)
        far_cutter = arrow_cutter.translate(far_rot)
        near_body = arrow_body.translate(near_rot)
        near_cutter = arrow_cutter.translate(near_rot)
        indicator = (
            indicator.union(far_body).cut(far_cutter).union(near_body).cut(near_cutter)
        )

        solid = foots.union(wall1).union(wall2).union(wall3).union(indicator)
        solid = solid.rotate((0, 0, 0), (0, -1, 0), sketch_angle_deg)

        solid = solid.cut(
            body.bounding_solid.translate((-GENERIC_GAP, 0, -GENERIC_GAP))
        )

        self.solid = solid
        self.angle_deg = sketch_angle_deg + 90


class FocusTower:
    def __init__(self):
        tower_diameter = 50
        tower_wall_thickness = 2

        x_offset = -tower_diameter / 2
        floor_thickness = 6
        floor_length = CEIL_OUTER_MOUNT_X + 5 - x_offset
        floor_width = (CEIL_OUTER_MOUNT_Y + 5) * 2

        foot_width = 15

        solid = (
            cq.Workplane("XY")
            .box(
                floor_length,
                floor_width,
                floor_thickness,
                centered=(False, True, False),
            )
            .translate((x_offset, 0, FRAME_CEIL_Z))
        )

        # 足
        foot = (
            cq.Workplane("XY")
            .box(
                60 + 10,
                foot_width,
                floor_thickness,
                centered=(False, True, False),
            )
            .translate((-10, 0, FRAME_CEIL_Z))
        )
        solid = solid.union(
            foot.rotate((0, 0, 0), (0, 0, 1), 45).translate(
                (CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y, 0)
            )
        )
        solid = solid.union(
            foot.rotate((0, 0, 0), (0, 0, 1), -45).translate(
                (CEIL_OUTER_MOUNT_X, -CEIL_OUTER_MOUNT_Y, 0)
            )
        )

        # 穴開け
        hole_points = [
            (CEIL_INNER_MOUNT_X, 0),
            (CEIL_OUTER_MOUNT_X, -CEIL_OUTER_MOUNT_Y),
            (CEIL_OUTER_MOUNT_X, CEIL_OUTER_MOUNT_Y),
        ]
        solid = (
            solid.faces(">Z")
            .workplane(origin=(0, 0, 0))
            .pushPoints(hole_points)
            .circle(M3_HOLE_DIAMETER / 2)
            .cutBlind(-30)
            .edges(">Z and %circle")
            .chamfer(2)
            .edges("<Z and %circle")
            .chamfer(CHAMFER)
            .edges("|Z and (not <<X)")
            .fillet(foot_width / 2 - 0.1)
            .edges("|Z and <<X")
            .fillet(floor_width / 2 - 0.1)
        )

        # タワーの壁
        r = tower_diameter / 2
        arrow_revolve_deg = 10
        wall_revolve_deg = 360 - 90 + arrow_revolve_deg
        verts = [
            (r, 0),
            (r, FRAME_CEIL_Z),
            (r - tower_wall_thickness, FRAME_CEIL_Z),
            (r - tower_wall_thickness, 0),
        ]
        solid = solid.union(
            cq.Workplane("YZ")
            .polyline(verts)
            .close()
            .revolve(wall_revolve_deg, (0, 0, 0), (0, 1, 0))
            .rotate((0, 0, 0), (0, 0, 1), -wall_revolve_deg / 2 + 90)
        )

        # 矢印
        arrow_r_space = 5
        arrow_root_height = (r - arrow_r_space) * 3 // 2
        arrow_edge_height = 1
        verts = [
            (arrow_r_space, 0),
            (arrow_r_space, arrow_edge_height),
            (r, arrow_root_height),
            (r, 0),
        ]
        arrow = (
            cq.Workplane("YZ")
            .polyline(verts)
            .close()
            .revolve(arrow_revolve_deg, (0, 0, 0), (0, 1, 0))
            .rotate((0, 0, 0), (0, 0, 1), -arrow_revolve_deg / 2 + 45)
        )
        for i in range(4):
            arrow_rotate_deg = i * (360 / 4)
            solid = solid.union(arrow.rotate((0, 0, 0), (0, 0, 1), arrow_rotate_deg))

        self.solid = solid


print("Generating mirror reflector segment...")
mirror_segment = MirrorSegment()

print("Generating mirror fastener...")
mirror_fastener = MirrorFastener()

print("Generating LED arm...")
led_arm = LedArm()

print("Generating focus indicator...")
focus_indicator = FocusIndicator(mirror_segment)
# focus_tower = FocusTower()

if True:
    preview_offset = 0
    mirror_reflector_solid = mirror_segment.reflector_solid.translate(
        (0, 0, preview_offset)
    )
    mirror_support_solid = mirror_segment.reflector_support.translate(
        (0, 0, preview_offset)
    )
    mirror_frame_solid = mirror_segment.frame_solid.translate(
        (0, 0, preview_offset * 2)
    )
    mirror_cap_solid = mirror_fastener.solid.translate((0, 0, preview_offset * 3))
    # led_arm_shaft_solid = led_arm.shaft2_solid.translate(
    #    (preview_offset, 0, preview_offset * 2)
    # )
    arm_joint1_solid = led_arm.joint1_solid.translate(
        (preview_offset, 0, preview_offset * 2)
    )
    arm_joint2_solid = led_arm.joint2_solid.translate(
        (preview_offset, 0, preview_offset * 2)
    )
    arm_left_link_solid = led_arm.left_link_solid.translate(
        (preview_offset, 0, preview_offset * 2)
    )
    arm_right_link_solid = led_arm.right_link_solid.translate(
        (preview_offset, 0, preview_offset * 2)
    )
    focus_indicator_solid = focus_indicator.solid
    # focus_tower_solid = focus_tower.solid
    show_object(mirror_reflector_solid, options={"color": "#eee"})
    # show_object(mirror_support_solid, options={"color": "#0f0"})
    show_object(mirror_frame_solid, options={"color": "#888"})
    show_object(mirror_cap_solid, options={"color": "#111"})
    show_object(arm_joint1_solid, options={"color": "#111"})
    # show_object(led_arm_shaft_solid, options={"color": "#111"})
    show_object(arm_joint2_solid, options={"color": "#111"})
    show_object(arm_left_link_solid, options={"color": "#888"})
    show_object(arm_right_link_solid, options={"color": "#888"})
    show_object(focus_indicator_solid, options={"color": "#84f"})
    # show_object(focus_tower_solid, options={"color": "#84f"})

    # 焦点の位置 (参考用)
    focus_marker = cq.Workplane("XY").box(2, 2, 2)
    show_object(focus_marker.translate(FOCUS_NEAR), options={"color": "#f00"})
    show_object(focus_marker.translate(FOCUS_FAR), options={"color": "#f00"})


def displace_reflector(solid):
    return solid.rotate((0, 0, 0), (0, -1, 0), -ELLIPSE_ANGLE).translate(
        (-ELLIPSE_LONG_RADIUS, 0, -REFLECTOR_FLOOR_Z)
    )


mirror_reflector_step = displace_reflector(mirror_segment.reflector_solid)
mirror_support_step = displace_reflector(mirror_segment.reflector_support)
if REFLECTOR_ENABLE_SUPPORT:
    mirror_reflector_step = mirror_reflector_step.union(mirror_support_step)
mirror_reflector_step.export(f"{STEP_OUT_DIR}/mirror_reflector.step")

mirror_frame_step = mirror_segment.frame_solid.rotate((0, 0, 0), (0, -1, 0), 180)
mirror_frame_step.export(f"{STEP_OUT_DIR}/mirror_frame.step")

mirror_fastener_step = mirror_fastener.solid
mirror_fastener_step.export(f"{STEP_OUT_DIR}/mirror_fastener.step")

arm_joint1_step = led_arm.joint1_solid.rotate((0, 0, 0), (0, -1, 0), 90)
arm_joint1_step.export(f"{STEP_OUT_DIR}/arm_joint1.step")
arm_joint2_step = led_arm.joint2_solid
arm_joint2_step.export(f"{STEP_OUT_DIR}/arm_joint2.step")
arm_left_link_step = led_arm.left_link_solid.rotate((0, 0, 0), (1, 0, 0), -90)
arm_left_link_step.export(f"{STEP_OUT_DIR}/arm_left_link.step")
arm_right_link_step = led_arm.right_link_solid.rotate((0, 0, 0), (1, 0, 0), 90)
arm_right_link_step.export(f"{STEP_OUT_DIR}/arm_right_link.step")

focus_indicator_step = focus_indicator.solid.rotate(
    (0, 0, 0), (0, 1, 0), focus_indicator.angle_deg
)
focus_indicator_step.export(f"{STEP_OUT_DIR}/focus_indicator.step")

# focus_tower_step = focus_tower.solid.rotate((0, 0, 0), (0, -1, 0), 180)
# focus_tower_step.export(f"{STEP_OUT_DIR}/focus_tower.step")

if False:
    show_object(mirror_reflector_step, options={"color": "#eee"})
    if REFLECTOR_ENABLE_SUPPORT:
        show_object(mirror_support_step, options={"color": "#0f0"})
