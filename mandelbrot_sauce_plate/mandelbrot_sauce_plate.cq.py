import cadquery as cq
import numpy as np

# 皿の床面の半径
FLOOR_R = 30

# 皿の底面から彫り込みの底までの厚み
BOTTOM_H = 3

# 彫り込みの深さ
IMAGE_DEPTH = 5

# 皿の底の高さ
FLOOR_H = BOTTOM_H + IMAGE_DEPTH

# マンデルブロ集合のピクセルサイズの目標値
PITCH_GOAL = 0.4

# マンデルブロ集合描画領域の半径
# (底の半径より少し狭くする)
IMAGE_R = FLOOR_R - PITCH_GOAL * 2

# マンデルブロ集合の解像度
NUM_PIXELS = int(np.ceil(IMAGE_R * 2 / PITCH_GOAL))

# マンデルブロ集合のピクセルサイズの実際値
PITCH = IMAGE_R * 2 / NUM_PIXELS

# 皿のフチの高さ
WALL_H = 5

# 皿のフチの水平方向の厚み
WALL_T = 5

# マンデルブロ集合の注目座標
RENDER_REAL = -0.5
RENDER_IMAG = 0

# マンデルブロ集合の描画半径
RENDER_R = 1.5

# マンデルブロ集合関数
def mandelbrot(a, b, max_iter):
    z = 0
    c = complex(a, b)
    for n in range(max_iter):
        if abs(z) > 2:
            return n
        z = z * z + c
    return max_iter

# 彫刻用の cutter 生成
def build(ix0, ix1):
    # 描画範囲の算出
    real0 = RENDER_REAL - RENDER_R
    real1 = RENDER_REAL + RENDER_R
    imag0 = RENDER_IMAG - RENDER_R
    imag1 = RENDER_IMAG + RENDER_R

    x0 = -IMAGE_R + IMAGE_R * 2 * ix0 / NUM_PIXELS
    x1 = -IMAGE_R + IMAGE_R * 2 * ix1 / NUM_PIXELS
    x_thick = x1 - x0
    x_middle = (x0 + x1) / 2
    a = RENDER_REAL + RENDER_R * x_middle / IMAGE_R

    if ix0 + 1 == ix1:
        # 縦一列分の描画
        
        # 縦方向の描画範囲の算出
        # (丸い皿の底から飛び出さないように)
        r = IMAGE_R - PITCH / 2
        iy0 = int(np.floor(-np.sqrt(r**2 - x_middle**2) / PITCH)) - 1
        iy1 = -iy0
        iy0 += int(NUM_PIXELS // 2)
        iy1 += int(NUM_PIXELS // 2)
        y0 = -IMAGE_R + IMAGE_R * 2 * iy0 / NUM_PIXELS
        y1 = -IMAGE_R + IMAGE_R * 2 * iy1 / NUM_PIXELS

        # 断面座標の生成
        max_iter = 100
        z_last = -1
        verts = [(y0, 5)]
        for iy in range(iy0, iy1):
            y = -IMAGE_R + IMAGE_R * 2 * iy / NUM_PIXELS
            b = RENDER_IMAG + RENDER_R * (y + PITCH / 2) / IMAGE_R
            z = (mandelbrot(a, b, max_iter)) / max_iter
            z = 1 - (1 - z) ** 4
            z = -IMAGE_DEPTH * z
            if z != z_last:
                verts.append((y, z_last))
                verts.append((y, z))
            z_last = z
        verts.append((y1, z_last))
        verts.append((y1, 5))
        
        # 押し出しで断面を生成
        return (
            cq.Workplane("YZ")
            .polyline(verts)
            .close()
            .extrude(x_thick)
            .translate((x0, 0, 0))
        )
    else:
        # 分割-->統合
        ixc = (ix0 + ix1) // 2
        sub0 = build(ix0, ixc)
        sub1 = build(ixc, ix1)
        return sub0.union(sub1)

# 皿の生成
dish_verts = [
    (0, 0),
    (0, FLOOR_H),
    (FLOOR_R, FLOOR_H),
    (FLOOR_R + WALL_H * 1.5, FLOOR_H + WALL_H),
    (FLOOR_R + WALL_H * 1.5 + WALL_T, FLOOR_H + WALL_H),
    (FLOOR_R + WALL_H * 0.5 - FLOOR_H, 0),
]
dish = (
    cq.Workplane("YZ")
    .polyline(dish_verts)
    .close()
    .revolve(360, (0, 0, 0), (0, 1, 0))
    .faces(">Z")
    .fillet(1)
    .faces("<Z")
    .fillet(5)
)

# マンデルブロ集合を彫り込むための cutter を生成
cutter = build(0, NUM_PIXELS).translate((0, 0, FLOOR_H + 0.4))

#show_object(cutter)

# 彫り込み
dish = dish.cut(cutter)

# 出力
show_object(dish)
dish.export("mandelbrot_sauce_plate.step")
