"""Download missing pinned source files; existing files must pass SHA256 verification."""
from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
 for r in json.loads((ROOT/'source_manifest.json').read_text()):
  if r['status']!='downloaded_valid':continue
  p=ROOT/r['local_path'];p.parent.mkdir(parents=True,exist_ok=True)
  if not p.exists():
   with urllib.request.urlopen(r['url']) as f:p.write_bytes(f.read())
  assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],f'Release bytes changed: {p}; investigate, do not silently update pin.'
  print('verified',p.name)
