from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: extract_scmarkeragent_zip.py <zip_path> <out_dir>")
        return 2

    zip_path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    report = []
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            item = {
                "filename": info.filename,
                "file_size": info.file_size,
                "compress_size": info.compress_size,
                "status": "pending",
            }
            try:
                if info.is_dir():
                    (out_dir / info.filename).mkdir(parents=True, exist_ok=True)
                    item["status"] = "dir"
                else:
                    target = out_dir / info.filename
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(info, "r") as src, open(target, "wb") as dst:
                        while True:
                            chunk = src.read(1024 * 1024)
                            if not chunk:
                                break
                            dst.write(chunk)
                    item["status"] = "ok"
            except Exception as exc:  # noqa: BLE001
                item["status"] = "error"
                item["error"] = repr(exc)
            report.append(item)

    report_path = out_dir / "extract_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=True, indent=2), encoding="utf-8")

    errors = [x for x in report if x["status"] == "error"]
    print(f"entries={len(report)} errors={len(errors)} report={report_path}")
    if errors:
        for err in errors[:10]:
            print(f"ERROR {err['filename']}: {err['error']}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
