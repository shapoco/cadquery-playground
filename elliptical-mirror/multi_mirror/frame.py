import cadquery as cq
import math

# 嵌合用の隙間 [mm]
GENERIC_GAP = 0.5

# 面取りのサイズ [mm]
CHAMFER = 0.5

# 固定用キャップの穴の間隔 [mm]
CAP_MOUNT_DISTANCE_X = 50
CAP_MOUNT_DISTANCE_Y = 30

M3_TAP_HOLE_DIAMETER = 2.5
M3_HOLE_DIAMETER = 3 + GENERIC_GAP

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

bridge_moved = bridge_main_solid.translate((0, CAP_MOUNT_DISTANCE_X / 2 - 5, 0))
show_object(bridge_moved, options={"color": "#888"})
show_object(bridge_moved.mirror("XZ"), name="mirror", options={"color": "#888"})

for p in sub_beam_hole_points:
    show_object(bridge_sub_solid.translate((p[0], 0, p[1])), options={"color": "#111"})


bridge_main_step = bridge_main_solid.rotate((0, 0, 0), (1, 0, 0), 90)
bridge_main_step.export(f"{STEP_OUT_DIR}/bridge_main.step")

bridge_sub_step = bridge_sub_solid.rotate((0, 0, 0), (1, 0, 0), 90)
bridge_sub_step.export(f"{STEP_OUT_DIR}/bridge_sub.step")
