#!/usr/bin/env python3
"""Check frozen source hashes, CLI coverage, source links and Python syntax.

This reads text/metadata only and does not compile or execute a model. An edited
source tree intentionally fails the frozen-identity check; it is a new version.
"""
from pathlib import Path
import ast,hashlib,json,re

def main():
    """Verify the published frozen-source contract, reporting precise failures."""
    r=Path(__file__).resolve().parents[1];m=json.loads((r/'validation/source_identity.json').read_text())
    for name,h in m['numerical_SHA256'].items():
        assert hashlib.sha256((r/name).read_bytes()).hexdigest()==h,('source modified',name)
    expected=hashlib.sha256(''.join(m['numerical_SHA256'][k] for k in sorted(m['numerical_SHA256'])).encode()).hexdigest()
    assert expected==m['source_ID']
    actual=set(re.findall(r'"(--[A-Za-z0-9-]+)"',(r/'src/main.c').read_text()))
    documented={x['option'] for x in json.loads((r/'docs/parameters.json').read_text())}
    assert actual==documented,('CLI documentation mismatch',actual^documented)
    for x in json.loads((r/'docs/function_contracts.json').read_text()):
        lines=(r/x['file']).read_text().splitlines();line=lines[x['line']-1]
        assert re.search(r'\b'+re.escape(x['name'])+r'\s*\(',line),('source link mismatch',x)
    for p in r.rglob('*.py'):
        if not any(z.startswith('build') for z in p.relative_to(r).parts):ast.parse(p.read_text(),filename=str(p))
    missing=[]
    for p in [r/'README.md',r/'THIRD_PARTY.md',*(r/'docs').glob('*.md'),r/'methodology/README.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            path=link.split('#')[0]
            if path and not (p.parent/path).exists():missing.append((str(p.relative_to(r)),link))
    assert not missing,('missing documentation targets',missing)
    frozen=r/'validation/release_files_SHA256.json'
    if frozen.exists():
        for name,h in json.loads(frozen.read_text())['files_SHA256'].items():
            assert hashlib.sha256((r/name).read_bytes()).hexdigest()==h,('release file modified',name)
    assert hashlib.sha256((r/'data/igrf14coeffs.txt').read_bytes()).hexdigest()=='8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c'
    print('RELEASE_TEXT_CHECK_PASS',m['source_ID'],'CLI_options='+str(len(actual)))

if __name__=='__main__':main()
