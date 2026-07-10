"""Convert TLS feature PDFs to JPEG and stitch into 4x5 JPG grids"""
import fitz, os, glob, sys
from PIL import Image

BASE = sys.argv[1] if len(sys.argv) > 1 else r"E:\GBM\results\tls_official_cut01"
FEATURES = [
    "tls_score_official",
    "tls_region_official",
    "plasma_b_cells",
    "t_cells",
    "bt_codistribution",
    "lc50_sig",
    "spatial_niche",
]
COLS, ROWS = 4, 5
N_PER_GRID = COLS * ROWS
TW, TH = 600, 600

for feat in FEATURES:
    print(f"\n=== {feat} ===")
    pdfs = sorted(glob.glob(f"{BASE}/*/{feat}.pdf"))
    print(f"Found {len(pdfs)} PDFs")

    # Convert PDF -> JPEG (in memory, no intermediate files)
    imgs = []
    for pdf_path in pdfs:
        doc = fitz.open(pdf_path)
        pix = doc[0].get_pixmap(dpi=200)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        doc.close()
        imgs.append(img)
    print(f"Converted {len(imgs)} PDFs")

    # Stitch into 4x5 JPG grids
    for g in range(0, len(imgs), N_PER_GRID):
        batch = imgs[g:g + N_PER_GRID]
        gi = g // N_PER_GRID + 1

        cells = []
        for img in batch:
            img.thumbnail((TW, TH), Image.LANCZOS)
            canvas = Image.new("RGB", (TW, TH), "white")
            canvas.paste(img, ((TW - img.width) // 2, (TH - img.height) // 2))
            cells.append(canvas)
        while len(cells) < N_PER_GRID:
            cells.append(Image.new("RGB", (TW, TH), "white"))

        grid = Image.new("RGB", (COLS * TW, ROWS * TH), "white")
        for i, img in enumerate(cells):
            grid.paste(img, ((i % COLS) * TW, (i // COLS) * TH))

        out = f"{BASE}/{feat}_grid_{gi:02d}.jpg"
        grid.save(out, quality=85)
        print(f"  Grid {gi}: {len(batch)} imgs -> {out}")

print("\nDone!")
