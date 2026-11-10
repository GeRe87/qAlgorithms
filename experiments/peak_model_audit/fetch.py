"""Fetch only the fixed public audit inputs and verify their published checksums."""
import hashlib
import json
from pathlib import Path
import sys
import urllib.request

DATASETS = [
    ("injection1", "20171016_POOL_POS_1_105-134.mzML", "4130309a9283c8b492660dbad1fd3670"),
    ("injection19", "20171016_POOL_POS_3_105-134.mzML", "8bce27a1fdbe2ea0f0659d398b0bfaaf"),
]

if __name__ == "__main__":
    destination = Path(sys.argv[1])
    destination.mkdir(parents=True, exist_ok=True)
    manifest = []
    for label, filename, expected in DATASETS:
        url = "https://zenodo.org/records/20729183/files/" + filename + "?download=1"
        target = destination / (label + ".mzML")
        with urllib.request.urlopen(url, timeout=120) as response:
            payload = response.read()
        actual = hashlib.md5(payload).hexdigest()
        if actual != expected:
            raise RuntimeError(f"Published checksum mismatch for {filename}")
        target.write_bytes(payload)
        manifest.append(dict(label=label, filename=filename, url=url, md5=actual,
                             sha256=hashlib.sha256(payload).hexdigest(), bytes=len(payload)))
    (destination / "input-manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
