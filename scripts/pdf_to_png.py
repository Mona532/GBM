"""Convert all tls_score_official.pdf to PNG"""
import fitz  # PyMuPDF
import os, glob

BASE = "E:/GBM/results/tls_official_relaxed"
pdfs = sorted(glob.glob(f"{BASE}/*/tls_score_official.pdf"))
print(f"Found {len(pdfs)} PDFs")

for pdf_path in pdfs:
    sample = os.path.basename(os.path.dirname(pdf_path))
    png_path = os.path.join(os.path.dirname(pdf_path), "tls_score_official.png")
    try:
        doc = fitz.open(pdf_path)
        page = doc[0]
        pix = page.get_pixmap(dpi=200)
        pix.save(png_path)
        doc.close()
        print(f"OK: {sample}")
    except Exception as e:
        print(f"ERR: {sample} - {e}")

print("Done")
