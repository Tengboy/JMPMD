"""Check a freshly extracted release (Python >=3.9, standard library only)."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]


def read_table(name):
    with (ROOT / 'results/published' / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def main():
    sums = ROOT / 'SHA256SUMS.txt'
    checked = 0
    if sums.exists():
        for line in sums.read_text('utf-8').splitlines():
            expected, rel = line.split('  ', 1)
            p = ROOT / rel
            assert p.is_file(), f'Missing file: {rel}'
            assert hashlib.sha256(p.read_bytes()).hexdigest() == expected, f'Changed file: {rel}'
            checked += 1
    else:
        raise AssertionError('SHA256SUMS.txt is missing; this is not a finalized release.')
    for p in ROOT.rglob('*'):
        if not p.is_file():
            continue
        assert p.suffix.lower() not in {'.mat', '.fig', '.log', '.asv', '.zip'}, f'Generated/binary trial file: {p}'
        if p.suffix in {'.m', '.json', '.csv', '.md', '.py', '.tex'}:
            text = p.read_text('utf-8-sig')
            assert not re.search(r'[A-Za-z]:[/\\]', text), f'Machine-specific drive path: {p}'
    summary = read_table('summary.csv')
    conv = read_table('convergence.csv')
    assert len(summary) == 232 and len(conv) == 800
    assert all(r['Strict'] == '1' and r['AngleSolver'] == 'local' and r['Pending'] == '0' for r in summary)
    assert {r['Group'] for r in summary} == {'A1', 'A2', 'A3', 'A4', 'B1', 'B3'}
    keys = [(r['Group'], r['Kind'], r['Series'], float(r['X'])) for r in summary]
    assert len(keys) == len(set(keys))
    spec = json.loads((ROOT / 'metadata/figure_layout.json').read_text('utf-8'))
    assert len(spec) == 7
    curves = points = 0
    for fig in spec:
        for panel in fig['panels']:
            for c in panel['curves']:
                assert 'y' not in c, 'Plot layout must read data from CSV, not store hidden substitute values.'
                if c['table'] == 'summary':
                    rows = [r for r in summary if r['Kind'] == c['kind'] and r['Group'] == c['group'] and r['Series'] == c['series']]
                    xfield = 'X'
                else:
                    rows = [r for r in conv if r['Kind'] == c['kind'] and float(r['SNR']) == c['snr']]
                    xfield = 'Iteration'
                for x in c['expected_x']:
                    chosen = [r for r in rows if float(r[xfield]) == x]
                    assert len(chosen) == 1, (fig['name'], c['label'], x)
                    y = float(chosen[0][c['metric']])
                    assert math.isfinite(y) and (not panel['log'] or y > 0)
                    points += 1
                curves += 1
        for ext in ('pdf', 'png'):
            assert (ROOT / 'results/published/figures' / f'{fig["name"]}.{ext}').is_file()
    print(json.dumps({'passed': True, 'hashed_files': checked, 'main_figures': 7,
                      'curves': curves, 'points': points, 'summary_rows': len(summary),
                      'convergence_rows': len(conv), 'raw_mat_files': 0}, indent=2))


if __name__ == '__main__':
    main()
