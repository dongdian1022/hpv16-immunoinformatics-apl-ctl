from PIL import Image, ImageDraw, ImageFont
import os

# ══════════════════════════════════════════════════════════════════
# 输入文件（修改为你的实际路径）
# ══════════════════════════════════════════════════════════════════
PATH_A = "Figure3A.jpg"
PATH_B = "Figure3B.jpg"
PATH_C = "Figure3C.jpg"     # 已处理好的Figure3C（含标注和比例尺）

OUTPUT_PNG  = "Figure3.png"
OUTPUT_TIFF = "Figure3.tiff"
OUTPUT_PDF  = "Figure3.pdf"
OUTPUT_JPG  = "Figure3.jpg"

# ══════════════════════════════════════════════════════════════════
# 布局设计
# 上排：A（折线图）+ B（柱状图），左右并排
# 下排：C（显微图），居中铺满
#
# 目标输出宽度：5400px（18 inches @ 300 dpi，期刊双栏宽度）
# ══════════════════════════════════════════════════════════════════
TARGET_W   = 5400
GAP        = 60     # panel 间距（px）
MARGIN     = 80     # 图四周边距
TAG_SIZE   = 130    # A/B/C panel tag 字号
TAG_OFFSET = 20     # tag 距左上角的偏移

# ── 加载图片 ────────────────────────────────────────────────────
img_a = Image.open(PATH_A).convert('RGB')
img_b = Image.open(PATH_B).convert('RGB')
img_c = Image.open(PATH_C).convert('RGB')

# ── 上排宽度分配：A:B = 原始宽度比例，高度统一 ─────────────────
top_area_w  = TARGET_W - MARGIN * 2 - GAP
ratio_ab    = img_a.width / img_b.width      # ≈ 1.18
w_a_top     = int(top_area_w * ratio_ab / (1 + ratio_ab))
w_b_top     = top_area_w - w_a_top

# A 和 B 缩放到相同高度（以 A 的高宽比为基准）
h_top       = int(w_a_top * img_a.height / img_a.width)
h_b_scaled  = int(w_b_top * img_b.height / img_b.width)
# 取两者较小值保证对齐
h_top       = min(h_top, h_b_scaled)

img_a_r = img_a.resize((w_a_top, h_top), Image.LANCZOS)
img_b_r = img_b.resize((int(h_top * img_b.width / img_b.height), h_top), Image.LANCZOS)
w_b_top = img_b_r.width   # 实际宽度（略有浮点修正）

# ── 下排：C 等比缩放到全宽 ────────────────────────────────────
w_c         = TARGET_W - MARGIN * 2
h_c         = int(w_c * img_c.height / img_c.width)
img_c_r     = img_c.resize((w_c, h_c), Image.LANCZOS)

# ── 计算总画布高度 ──────────────────────────────────────────────
total_h = MARGIN + h_top + GAP + h_c + MARGIN

# ── 创建白色画布 ────────────────────────────────────────────────
canvas = Image.new('RGB', (TARGET_W, total_h), (255, 255, 255))

# ── 粘贴各 panel ────────────────────────────────────────────────
x_a = MARGIN
x_b = MARGIN + w_a_top + GAP
y_top = MARGIN
y_c   = MARGIN + h_top + GAP

canvas.paste(img_a_r, (x_a, y_top))
canvas.paste(img_b_r, (x_b, y_top))
canvas.paste(img_c_r, (MARGIN, y_c))

# ── 加 panel tag（A / B / C），用 textbbox 精确定位 ─────────────
draw = ImageDraw.Draw(canvas)

try:
    font_tag = ImageFont.truetype(
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc", TAG_SIZE)
except:
    font_tag = ImageFont.load_default()

def add_tag(draw, text, x, y, font):
    """在 (x,y) 处添加加粗 panel tag，白色背景防遮挡"""
    bb   = draw.textbbox((0, 0), text, font=font)
    tw   = bb[2] - bb[0]
    th   = bb[3] - bb[1]
    pad  = 8
    # 白色背景框
    draw.rectangle([x-pad, y-pad, x+tw+pad, y+th+pad],
                   fill=(255, 255, 255))
    draw.text((x, y), text, fill=(0, 0, 0), font=font)

add_tag(draw, "A", x_a + TAG_OFFSET,       y_top + TAG_OFFSET, font_tag)
add_tag(draw, "B", x_b + TAG_OFFSET,       y_top + TAG_OFFSET, font_tag)
add_tag(draw, "C", MARGIN + TAG_OFFSET,    y_c   + TAG_OFFSET, font_tag)

print(f"画布尺寸: {canvas.size}")
print(f"  A panel: ({x_a}, {y_top})  {img_a_r.size}")
print(f"  B panel: ({x_b}, {y_top})  {img_b_r.size}")
print(f"  C panel: ({MARGIN}, {y_c})  {img_c_r.size}")

# ══════════════════════════════════════════════════════════════════
# 保存（均 300 dpi）
# ══════════════════════════════════════════════════════════════════
canvas.save(OUTPUT_PNG,  dpi=(300, 300))
print(f"✅ PNG   {os.path.getsize(OUTPUT_PNG)/1024/1024:.1f} MB")

canvas.save(OUTPUT_TIFF, dpi=(300, 300), compression='tiff_lzw')
print(f"✅ TIFF  {os.path.getsize(OUTPUT_TIFF)/1024/1024:.1f} MB")

canvas.save(OUTPUT_JPG,  dpi=(300, 300), quality=95, subsampling=0)
print(f"✅ JPG   {os.path.getsize(OUTPUT_JPG)/1024/1024:.1f} MB")

canvas.save(OUTPUT_PDF,  resolution=300)
print(f"✅ PDF   {os.path.getsize(OUTPUT_PDF)/1024/1024:.1f} MB")

print("完成！")

# pip install pillow
