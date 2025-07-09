import cadquery as cq
import math

# 嵌合用の隙間 [mm]
GENERIC_GAP = 0.5
SHAFT_GAP = 0.2

# 面取りのサイズ [mm]
CHAMFER = 0.5

# 固定用キャップの穴の間隔 [mm]
CAP_MOUNT_DISTANCE_X = 50
CAP_MOUNT_DISTANCE_Y = 30

M3_TAP_HOLE_DIAMETER = 2.5
M3_HOLE_DIAMETER = 3 + GENERIC_GAP

M6_NUT_DIAMETER = 10
M6_NUT_THICKNESS = 5

BRIDGE_HEIGHT = 40

BEAM_HEIGHT = 10
BEAM_WIDTH = 5
BEAM_LENGTH = 200

PILLAR_THICKNESS = 5

HOLE_STAGE_THICKNESS = 3
HOLE_STAGE_WIDTH = 8 + HOLE_STAGE_THICKNESS * 2
HOLE_STAGE_LENGTH = 10

STEP_OUT_DIR = "./step"


bridge_main_solid = (
    cq.Workplane("XY")
    .box(BEAM_LENGTH, BEAM_WIDTH, BEAM_HEIGHT, centered=(True, False, False))
    .translate((0, -BEAM_WIDTH, 0))
)
verts = [
    (0, 0),
    (10, 0),
    (5, 5),
    (0, 5),
]
joint = (
    cq.Workplane("XY")
    .polyline(verts)
    .close()
    .extrude(BEAM_HEIGHT)
    .translate((-BEAM_LENGTH / 2, 0, 0))
)
bridge_main_solid = bridge_main_solid.union(joint)
bridge_main_solid = bridge_main_solid.union(joint.mirror("YZ"))
bridge_main_solid = (
    bridge_main_solid.faces(">X")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(0, 5)])
    .circle(M3_TAP_HOLE_DIAMETER / 2)
    .cutBlind(-5)
    .faces("<X")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(0, 5)])
    .circle(M3_TAP_HOLE_DIAMETER / 2)
    .cutBlind(-7.5)
)

bridge_main_solid = bridge_main_solid.union(
    bridge_main_solid.translate((0, 0, BRIDGE_HEIGHT - BEAM_HEIGHT))
)

pillar = (
    cq.Workplane("XY")
    .box(
        PILLAR_THICKNESS,
        10,
        BRIDGE_HEIGHT - BEAM_HEIGHT * 2,
        centered=(True, False, False),
    )
    .translate((BEAM_LENGTH / 2 - PILLAR_THICKNESS / 2, -BEAM_WIDTH, BEAM_HEIGHT))
)
bridge_main_solid = bridge_main_solid.union(pillar)
bridge_main_solid = bridge_main_solid.union(pillar.mirror("YZ"))

NUM_SUJIKAI = 5
SUJIKAI_THICKNESS = 10
SUJIKAI_WIDTH = (BEAM_LENGTH - SUJIKAI_THICKNESS) / NUM_SUJIKAI

x0 = 0
x1 = SUJIKAI_WIDTH / 2
y0 = BEAM_HEIGHT
y1 = BRIDGE_HEIGHT - BEAM_HEIGHT
verts = [
    (x0 - SUJIKAI_THICKNESS / 2, y1),
    (x0 + SUJIKAI_THICKNESS / 2, y1),
    (x1 + SUJIKAI_THICKNESS / 2, y0),
    (x1 - SUJIKAI_THICKNESS / 2, y0),
]
sujikai = cq.Workplane("XZ").polyline(verts).close().extrude(BEAM_WIDTH)
sujikai = sujikai.union(sujikai.mirror("YZ"))

x_offset = -BEAM_LENGTH / 2 + SUJIKAI_WIDTH / 2 + SUJIKAI_THICKNESS / 2
for i in range(NUM_SUJIKAI):
    bridge_main_solid = bridge_main_solid.union(
        sujikai.translate((x_offset + i * SUJIKAI_WIDTH, 0, 0))
    )

sub_beam_hole_points = [
    (-SUJIKAI_WIDTH * 3 / 2, BEAM_HEIGHT / 2),
    (SUJIKAI_WIDTH * 3 / 2, BEAM_HEIGHT / 2),
    (0, BRIDGE_HEIGHT - BEAM_HEIGHT / 2),
]
bridge_main_solid = (
    bridge_main_solid.faces("<Y")
    .workplane(origin=(0, 0, 0))
    .pushPoints(sub_beam_hole_points)
    .circle(M3_HOLE_DIAMETER / 2)
    .cutBlind(-100)
    .faces("-Y")
    .edges("%circle")
    .chamfer(CHAMFER)
    .faces("+Y")
    .edges("%circle")
    .chamfer(2)
    .faces(">>X or <<X")
    .edges("%circle")
    .chamfer(CHAMFER)
    # .edges("%circle")
    # .edges("<<Y[0]")
    # .chamfer(2)
    .edges("%line and >>Y and (>>Z or <<Z)")
    .chamfer(1)
)

hole_stage = (
    cq.Workplane("XY")
    .box(
        HOLE_STAGE_WIDTH,
        HOLE_STAGE_LENGTH,
        HOLE_STAGE_THICKNESS,
        centered=(True, False, False),
    )
    .faces(">Z")
    .workplane()
    .pushPoints([(0, 5)])
    .circle(M3_HOLE_DIAMETER / 2)
    .cutBlind(-100)
    .edges("%circle")
    .chamfer(CHAMFER)
    .translate((0, 0, BEAM_HEIGHT - HOLE_STAGE_THICKNESS))
)

verts = [
    (0, 0),
    (0, BEAM_HEIGHT - HOLE_STAGE_THICKNESS),
    (HOLE_STAGE_LENGTH, BEAM_HEIGHT - HOLE_STAGE_THICKNESS),
]

hole_stage_foot = (
    cq.Workplane("YZ")
    .polyline(verts)
    .close()
    .extrude(HOLE_STAGE_THICKNESS)
    .translate((HOLE_STAGE_WIDTH / 2 - HOLE_STAGE_THICKNESS, 0, 0))
)
hole_stage = hole_stage.union(hole_stage_foot)
hole_stage = hole_stage.union(hole_stage_foot.mirror("YZ"))

hole_stage = (
    hole_stage.faces(">Y").edges("|Z").chamfer(2).edges("|Y and (not >>Z)").chamfer(1)
)

bridge_main_solid = bridge_main_solid.union(
    hole_stage.translate((CAP_MOUNT_DISTANCE_Y / 2, 0, 0))
)
bridge_main_solid = bridge_main_solid.union(
    hole_stage.translate((-CAP_MOUNT_DISTANCE_Y / 2, 0, 0))
)


y1 = (CAP_MOUNT_DISTANCE_X - (5 * 2) - (BEAM_WIDTH * 2)) / 2
verts = [
    (0, 0),
    (2.5, 0),
    (2.5, y1 - 10),
    (5, y1 - 5),
    (5, y1),
    (0, y1),
]
bridge_sub_solid = (
    cq.Workplane("XY")
    .polyline(verts)
    .close()
    .extrude(BEAM_HEIGHT)
    .translate((0, 0, -BEAM_HEIGHT / 2))
)
bridge_sub_solid = bridge_sub_solid.union(bridge_sub_solid.mirror("YZ"))
bridge_sub_solid = bridge_sub_solid.union(bridge_sub_solid.mirror("XZ"))
bridge_sub_solid = (
    bridge_sub_solid.faces(">Y")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(0, 0)])
    .circle(M3_TAP_HOLE_DIAMETER / 2)
    .cutBlind(-7.5)
    .faces("<Y")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(0, 0)])
    .circle(M3_TAP_HOLE_DIAMETER / 2)
    .cutBlind(-7.5)
    .edges("%circle")
    .chamfer(CHAMFER)
)

GUIDE_SHAFT_DIAMETER = 4
SCREW_DIAMETER = 6

TSUMAMI_DIAMETER = 20
TSUMAMI_HEIGHT = 20

SHAFT_GRIP_BASE_THICKNESS = 4
SHAFT_GRIP_SHAFT_OFFSET = 10
SCREW_GUIDE_WIDTH = 20
SCREW_GUIDE_THICKNESS = 12
SCREW_GUIDE_HEIGHT = 10

elevator_grip_solid = (
    cq.Workplane("YZ")
    .box(
        CAP_MOUNT_DISTANCE_Y + 20,
        BRIDGE_HEIGHT,
        SHAFT_GRIP_BASE_THICKNESS,
        centered=(True, False, False),
    )
    .faces(">X")
    .workplane(origin=(0, 0, 0))
    .rect(TSUMAMI_DIAMETER + 5, TSUMAMI_HEIGHT, centered=(True, False))
    .cutBlind(-100)
)
hole_points = [
    (CAP_MOUNT_DISTANCE_Y / 2 + 5, 5),
    (-CAP_MOUNT_DISTANCE_Y / 2 - 5, 5),
    (CAP_MOUNT_DISTANCE_Y / 2 + 5, BRIDGE_HEIGHT - 5),
    (-CAP_MOUNT_DISTANCE_Y / 2 - 5, BRIDGE_HEIGHT - 5),
]
elevator_grip_solid = (
    elevator_grip_solid.faces(">X")
    .workplane(origin=(0, 0, 0))
    .pushPoints(hole_points)
    .circle(M3_HOLE_DIAMETER / 2)
    .cutBlind(-100)
    .faces(">X")
    .edges("%circle")
    .chamfer(2)
    .faces("<X")
    .edges("%circle")
    .chamfer(CHAMFER)
    .edges("|X")
    .chamfer(1)
)

hole_x = SHAFT_GRIP_SHAFT_OFFSET - SHAFT_GRIP_BASE_THICKNESS
screw_guide = cq.Workplane("XY").box(
    hole_x,
    SCREW_GUIDE_WIDTH,
    SCREW_GUIDE_HEIGHT,
    centered=(False, True, False),
)
screw_guide = screw_guide.union(
    cq.Workplane("XY")
    .cylinder(SCREW_GUIDE_HEIGHT, SCREW_GUIDE_WIDTH / 2, centered=(True, True, False))
    .translate((hole_x, 0, 0))
)
screw_guide = screw_guide.cut(
    cq.Workplane("XY")
    .box(99, 99, 99, centered=(False, True, True))
    .translate((SCREW_GUIDE_THICKNESS, 0, 0))
)
screw_guide = (
    screw_guide.faces("<Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(hole_x, 0)])
    .circle((SCREW_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-100)
    .edges("%circle or >>X")
    .fillet(CHAMFER)
    .translate((SHAFT_GRIP_BASE_THICKNESS, 0, TSUMAMI_HEIGHT))
)

elevator_grip_solid = elevator_grip_solid.union(screw_guide)

hole_x = SHAFT_GRIP_SHAFT_OFFSET - SHAFT_GRIP_BASE_THICKNESS
shaft_guide = cq.Workplane("XY").box(
    hole_x,
    10,
    BRIDGE_HEIGHT - 20,
    centered=(False, True, False),
)
shaft_guide = shaft_guide.union(
    cq.Workplane("XY")
    .cylinder(BRIDGE_HEIGHT - 20, 5, centered=(True, True, False))
    .translate((hole_x, 0, 0))
)
shaft_guide = (
    shaft_guide.faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(hole_x, 0)])
    .circle((GUIDE_SHAFT_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-100)
    .edges("%circle")
    .chamfer(CHAMFER)
)
# cut_angle_deg = 90
# pole_guide = pole_guide.cut(
#     cq.Workplane("XZ")
#     .rect(GUIDE_POLE_DIAMETER * 2, 100, centered=(False, True))
#     .revolve(cut_angle_deg, (0, 0, 0), (0, 1, 0))
#     .rotate((0, 0, 0), (0, 0, 1), -cut_angle_deg / 2)
#     .translate((ELEVATOR_BASE_THICKNESS, 0, 0))
# )
shaft_guide = shaft_guide.translate(
    (SHAFT_GRIP_BASE_THICKNESS, CAP_MOUNT_DISTANCE_Y / 2 + 5, 10)
)
elevator_grip_solid = elevator_grip_solid.union(shaft_guide)
elevator_grip_solid = elevator_grip_solid.union(shaft_guide.mirror("XZ"))

# --------

BASE_FLOOR_LENGTH = 25
BASE_FLOOR_THICKNESS = 5
base_solid = cq.Workplane("XY").box(
    BASE_FLOOR_LENGTH,
    CAP_MOUNT_DISTANCE_X,
    BASE_FLOOR_THICKNESS,
    centered=(False, True, False),
)
elevator_mount_hole_points = [
    (20, -15),
    (20, 15),
]
base_solid = (
    base_solid.faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints(elevator_mount_hole_points)
    .circle(M3_HOLE_DIAMETER / 2)
    .cutBlind(-100)
    .faces(">Z")
    .edges("%circle")
    .chamfer(2)
    .faces("<Z")
    .edges("%circle")
    .chamfer(CHAMFER)
    .edges("|Z")
    .chamfer(1)
)

shaft_guide = (
    cq.Workplane("XY")
    .cylinder(5, 5, centered=(True, True, False))
    .faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(0, 0)])
    .circle((GUIDE_SHAFT_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-100)
    .translate(
        (
            SHAFT_GRIP_SHAFT_OFFSET,
            CAP_MOUNT_DISTANCE_X / 2 - 5,
            BASE_FLOOR_THICKNESS,
        )
    )
)
base_solid = base_solid.union(shaft_guide)
base_solid = base_solid.union(shaft_guide.mirror("XZ"))
screw_guide = (
    cq.Workplane("XY")
    .cylinder(5, 10, centered=(True, True, False))
    .translate((SHAFT_GRIP_SHAFT_OFFSET, 0, BASE_FLOOR_THICKNESS))
)
base_solid = base_solid.union(screw_guide)
base_solid = (
    base_solid.faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(SHAFT_GRIP_SHAFT_OFFSET, 0)])
    .circle((SCREW_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-100)
    .faces("<Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(SHAFT_GRIP_SHAFT_OFFSET, 0)])
    .polygon(6, M6_NUT_DIAMETER * 2 / math.sqrt(3) + GENERIC_GAP)
    .cutBlind(-M6_NUT_THICKNESS - GENERIC_GAP)
    .edges(">>Z")
    .chamfer(CHAMFER)
)

# --------

cap_solid = cq.Workplane("XY").box(
    10, CAP_MOUNT_DISTANCE_X - 10, 10, centered=(True, True, False)
)
cap_solid = cap_solid.union(
    cq.Workplane("XY")
    .cylinder(10, 5, centered=(True, True, False))
    .translate((0, (CAP_MOUNT_DISTANCE_X - 10) / 2, 0))
)
cap_solid = cap_solid.union(
    cq.Workplane("XY")
    .cylinder(10, 5, centered=(True, True, False))
    .translate((0, -(CAP_MOUNT_DISTANCE_X - 10) / 2, 0))
)
cap_solid = (
    cap_solid.faces("<Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints(
        [
            (0, (CAP_MOUNT_DISTANCE_X - 10) / 2),
            (0, -(CAP_MOUNT_DISTANCE_X - 10) / 2),
        ]
    )
    .circle((GUIDE_SHAFT_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-5)
    .pushPoints([(0, 0)])
    .circle((SCREW_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-100)
)
verts = [
    (-15, 10),
    (15, 10),
    (10, 5),
    (-10, 5),
]
cap_solid = (
    cap_solid.faces(">X")
    .workplane(origin=(0, 0, 0))
    .polyline(verts)
    .close()
    .cutBlind(-99)
    .edges()
    .chamfer(CHAMFER)
)
cap_solid = cap_solid.translate((SHAFT_GRIP_SHAFT_OFFSET, 0, 0))

# --------

tsumami_num_verts = 48
verts = []
for i in range(tsumami_num_verts):
    r = TSUMAMI_DIAMETER / 2
    if i % 2 == 0:
        r -= 0.5
    angle = i * (360 / tsumami_num_verts)
    rad = math.radians(angle)
    x = r * math.cos(rad)
    y = r * math.sin(rad)
    verts.append((x, y))
tsumami_solid = (
    cq.Workplane("XY")
    .polyline(verts)
    .close()
    .extrude(TSUMAMI_HEIGHT - 1)
    .edges(">>Z or <<Z")
    .chamfer(CHAMFER)
)
tsumami_solid = tsumami_solid.union(
    cq.Workplane("XY")
    .cylinder(1, SCREW_DIAMETER / 2 + 5, centered=(True, True, False))
    .translate((0, 0, TSUMAMI_HEIGHT - 1))
)
tsumami_solid = (
    tsumami_solid.faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints([(0, 0)])
    .circle((SCREW_DIAMETER + SHAFT_GAP) / 2)
    .cutBlind(-100)
    .faces("<Z")
    .workplane(origin=(0, 0, 0))
    .polygon(6, M6_NUT_DIAMETER * 2 / math.sqrt(3) + GENERIC_GAP)
    .cutBlind(-M6_NUT_THICKNESS)
    .edges("%circle and >>Z")
    .chamfer(CHAMFER)
)
tsumami_solid = tsumami_solid.translate((SHAFT_GRIP_SHAFT_OFFSET, 0, 0))


# --------

TOWER_BEAM_THICKNESS = 5
TOWER_CEIL_THICKNESS = 10
TOWER_HEIGHT = 100
TOWER_WIDTH = CAP_MOUNT_DISTANCE_X

tower_solid = (
    cq.Workplane("XY")
    .box(
        BASE_FLOOR_LENGTH,
        TOWER_WIDTH,
        TOWER_CEIL_THICKNESS,
        centered=(False, True, False),
    )
    .faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints(elevator_mount_hole_points)
    .circle(M3_TAP_HOLE_DIAMETER / 2)
    .cutBlind(-7.5)
    .faces(">Z")
    .edges("%circle")
    .chamfer(CHAMFER)
    .translate((0, 0, TOWER_HEIGHT - TOWER_CEIL_THICKNESS))
)

pillar = (
    cq.Workplane("XY")
    .box(
        BASE_FLOOR_LENGTH,
        TOWER_BEAM_THICKNESS,
        TOWER_HEIGHT,
        centered=(False, False, False),
    )
    .translate((0, TOWER_WIDTH / 2 - TOWER_BEAM_THICKNESS, 0))
)
tower_solid = tower_solid.union(pillar)
tower_solid = tower_solid.union(pillar.mirror("XZ"))

tower_solid = tower_solid.edges(">Y or <Y or >Z").chamfer(CHAMFER)

sujikai_height = (TOWER_HEIGHT - TOWER_CEIL_THICKNESS - TOWER_BEAM_THICKNESS) / 2
sujikai_width = TOWER_WIDTH - TOWER_BEAM_THICKNESS * 2
sujikai_thickness = 5
verts = [
    (-sujikai_width / 2, -sujikai_thickness / 2),
    (-sujikai_width / 2, sujikai_thickness / 2),
    (sujikai_width / 2, sujikai_height + sujikai_thickness / 2),
    (sujikai_width / 2, sujikai_height - sujikai_thickness / 2),
]
sujikai = cq.Workplane("YZ").polyline(verts).close().extrude(BASE_FLOOR_LENGTH)
sujikai = sujikai.union(sujikai.mirror("XZ"))
sujikai = sujikai.union(sujikai.mirror("XY"))
sujikai = sujikai.translate((0, 0, TOWER_BEAM_THICKNESS + sujikai_height))
tower_solid = tower_solid.union(sujikai)

hole_points = [
    (5, -(TOWER_WIDTH / 2 + 10)),
    (5, (TOWER_WIDTH / 2 + 10)),
    (BASE_FLOOR_LENGTH - 5, -(TOWER_WIDTH / 2 + 10)),
    (BASE_FLOOR_LENGTH - 5, (TOWER_WIDTH / 2 + 10)),
]
floor = (
    cq.Workplane("XY")
    .box(
        BASE_FLOOR_LENGTH,
        TOWER_WIDTH + 15 * 2,
        TOWER_BEAM_THICKNESS,
        centered=(False, True, False),
    )
    .faces(">Z")
    .workplane(origin=(0, 0, 0))
    .pushPoints(hole_points)
    .circle(M3_HOLE_DIAMETER / 2)
    .cutBlind(-7.5)
    .faces(">Z")
    .edges("%circle")
    .chamfer(2)
    .faces("<Z")
    .edges("%circle")
    .chamfer(CHAMFER)
    .edges("|Z")
    .chamfer(2)
    .edges("%line and >Z")
    .chamfer(CHAMFER)
)
tower_solid = tower_solid.union(floor)

# --------

preview_offset = 5

bridge_moved = bridge_main_solid.translate((0, CAP_MOUNT_DISTANCE_X / 2 - 5, 0))
show_object(bridge_moved, options={"color": "#888"})
show_object(bridge_moved.mirror("XZ"), options={"color": "#888"})

for p in sub_beam_hole_points:
    show_object(bridge_sub_solid.translate((p[0], 0, p[1])), options={"color": "#111"})

elevator_grip_moved = elevator_grip_solid.translate(
    (BEAM_LENGTH / 2 + preview_offset, 0, 0)
)
show_object(elevator_grip_moved, options={"color": "#111"})
show_object(elevator_grip_moved.mirror("YZ"), options={"color": "#111"})

base_moved = base_solid.translate(
    (BEAM_LENGTH / 2 + preview_offset, 0, -preview_offset * 4)
)
show_object(base_moved, options={"color": "#111"})
show_object(base_moved.mirror("YZ"), options={"color": "#111"})

cap_moved = cap_solid.translate(
    (BEAM_LENGTH / 2 + preview_offset, 0, BRIDGE_HEIGHT - 10 + preview_offset)
)
show_object(cap_moved, options={"color": "#111"})
show_object(cap_moved.mirror("YZ"), options={"color": "#111"})

tsumami_moved = tsumami_solid.translate(
    (BEAM_LENGTH / 2 + preview_offset, 0, -preview_offset)
)
show_object(tsumami_moved, options={"color": "#888"})
show_object(tsumami_moved.mirror("YZ"), options={"color": "#888"})

tower_moved = tower_solid.translate(
    (BEAM_LENGTH / 2 + preview_offset, 0, -TOWER_HEIGHT - preview_offset * 5)
)
show_object(tower_moved, options={"color": "#888"})

# --------

bridge_main_step = bridge_main_solid.rotate((0, 0, 0), (1, 0, 0), 90)
bridge_main_step.export(f"{STEP_OUT_DIR}/bridge_main.step")

bridge_sub_step = bridge_sub_solid
bridge_sub_step.export(f"{STEP_OUT_DIR}/bridge_sub.step")

elevator_grip_step = elevator_grip_solid.rotate((0, 0, 0), (0, 1, 0), -90)
elevator_grip_step.export(f"{STEP_OUT_DIR}/elevator_grip.step")

base_step = base_solid
base_step.export(f"{STEP_OUT_DIR}/elevator_base.step")

cap_step = cap_solid
cap_step.export(f"{STEP_OUT_DIR}/elevator_cap.step")

tsumami_step = tsumami_solid
tsumami_step.export(f"{STEP_OUT_DIR}/elevator_tsumami.step")
