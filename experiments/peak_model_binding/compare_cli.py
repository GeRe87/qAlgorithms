"""Compare all synthetic mzML feature/group CSVs against separate baseline binary."""
import subprocess
import sys
from pathlib import Path
import tempfile
script=Path(__file__).resolve().parents[2]/'test/integration/test_mzml_cli.py'
with tempfile.TemporaryDirectory() as tmp:
    roots=[]
    for i,exe in enumerate(sys.argv[1:3]):
        root=Path(tmp)/str(i);root.mkdir()
        subprocess.run([sys.executable,str(script),str(Path(exe).resolve()),'--work-dir',str(root)],check=True)
        roots.append(next(root.iterdir()))
    files=[{str(p.relative_to(root)):p.read_bytes() for p in root.rglob('*.csv') if p.is_file()} for root in roots]
    assert files[0] and files[0]==files[1],'CLI output differential mismatch'
    print('Exact matching CLI CSVs:',len(files[0]))
