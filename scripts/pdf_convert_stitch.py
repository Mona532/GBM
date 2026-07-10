"""Convert tls_score_official.pdf to PNG and stitch into 4x5 grids"""
import fitz
import os, glob
from PIL import Image

BASE = r"E:\GBM\results\tls_official_relaxed"
pdfs = sorted(glob.glob(f"{BASE}/*/tls_score_official.pdf"))
print(f"Found {len(pdfs)} PDFs")

# Step 1: Convert to PNG
png_files = []
for i, pdf_path in enumerate(pdfs):
    sample_dir = os.path.dirname(pdf_path)
    sample = os.path.basename(sample_dir)
    png_path = os.path.join(sample_dir, "tls_score_official.png")
    if os.path.exists(png_path):
        png_files.append(png_path)
        continue
    try:
        doc = fitz.open(pdf_path)
        page = doc[0]
        pix = page.get_pixmap(dpi=200)
        pix.save(png_path)
        doc.close()
        png_files.append(png_path)
        if (i+1) % 20 == 0:
            print(f"  Converted {i+1}/{len(pdfs)}")
    except Exception as e:
        print(f"ERR {sample}: {e}")

print(f"Converted {len(png_files)} PNGs")

# Step 2: Stitch 4x5 grids
n_per_grid = 20
cols, rows = 4, 5
thumb_w, thumb_h = 600, 600

for g in range(0, len(png_files), n_per_grid):
    batch = png_files[g:g+n_per_grid]
    grid_idx = g // n_per_grid + 1
    print(f"Grid {grid_idx}: {len(batch)} images...")

    # Pad to 20 if needed
    cell_imgs = []
    for f in batch:
        img = Image.open(f)
        img.thumbnail((thumb_w, thumb_h), Image.LANCZOS)
        # Center on white canvas
        canvas = Image.new('RGB', (thumb_w, thumb_h), 'white')
        x = (thumb_w - img.width) // 2
        y = (thumb_h - img.height) // 2
        canvas.paste(img, (x, y))
        cell_imgs.append(canvas)

    while len(cell_imgs) < n_per_grid:
        cell_imgs.append(Image.new('RGB', (thumb_w, thumb_h), 'white'))

    # Build grid
    grid_w = cols * thumb_w
    grid_h = rows * thumb_h
    grid_img = Image.new('RGB', (grid_w, grid_h), 'white')

    for i, img in enumerate(cell_imgs):
        r, c = i // cols, i % cols
        grid_img.paste(img, (c * thumb_w, r * thumb_h))

    out_path = f"{BASE}/tls_score_grid_{grid_idx:02d}.png"
    grid_img.save(out_path)
    print(f"  -> {out_path}")

print("Done!")
