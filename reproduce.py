import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download(url, path, sha256):
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={'User-Agent': 'ADMET-identifiability-reproduction'})
    with urllib.request.urlopen(request, timeout=120) as response, path.open('wb') as stream:
        shutil.copyfileobj(response, stream)
    if digest(path) != sha256:
        raise ValueError(f'Input checksum mismatch: {path.name}')


def matches(actual, expected, tolerance):
    if isinstance(expected, dict):
        return isinstance(actual, dict) and actual.keys() == expected.keys() and all(
            matches(actual[k], v, tolerance) for k, v in expected.items())
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(
            matches(a, b, tolerance) for a, b in zip(actual, expected))
    if isinstance(expected, float):
        return isinstance(actual, (int, float)) and math.isclose(actual, expected, rel_tol=0, abs_tol=tolerance)
    return type(actual) is type(expected) and actual == expected


def main():
    parser = argparse.ArgumentParser(description="Reproduce the study's analyses and figures from saved scores and selections.")
    parser.add_argument('--workdir', type=Path, default=Path('work'), help='New directory for inputs and outputs')
    parser.add_argument('--inputs', type=Path, help='Use a local input archive instead of downloading it')
    args = parser.parse_args()
    sources = json.loads((ROOT / 'data_sources.json').read_text(encoding='utf-8'))
    expected = json.loads((ROOT / 'expected_results.json').read_text(encoding='utf-8'))
    if args.inputs is None and not sources['frozen_inputs']['url']:
        parser.error('No input download URL is configured; supply a local archive with --inputs.')
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    for folder in ['scripts', 'reporting/analysis']:
        shutil.copytree(ROOT / folder, work / folder)
    for folder in ['reporting/tables', 'reporting/figures', 'logs']:
        (work / folder).mkdir(parents=True, exist_ok=True)
    archive = work / 'frozen_inputs.zip'
    if args.inputs:
        shutil.copyfile(args.inputs, archive)
    else:
        download(sources['frozen_inputs']['url'], archive, sources['frozen_inputs']['sha256'])
    if digest(archive) != sources['frozen_inputs']['sha256']:
        raise ValueError('Input archive checksum mismatch')
    with zipfile.ZipFile(archive) as package:
        allowed = {entry['path'] for entry in sources['frozen_inputs']['files']} | {'README.txt'}
        if set(package.namelist()) != allowed:
            raise ValueError('Unexpected files in the input archive')
        package.extractall(work)
    for entry in sources['frozen_inputs']['files']:
        if digest(work / entry['path']) != entry['sha256']:
            raise ValueError(f'Input checksum mismatch: {entry["path"]}')
    for entry in sources['source_data']:
        download(entry['url'], work / entry['path'], entry['sha256'])
    commands = [
        'scripts/analyze_expansionrx_v5_identifiability_lead_v1.py',
        'scripts/analyze_expansionrx_v6_strengthening_v1.py',
        'scripts/permutation_refinement_v6.py',
        'reporting/analysis/recompute_revision.py',
    ]
    plots = [
        'scripts/plot_selection_identification_v1.py',
        'reporting/analysis/plot_revision.py',
        'reporting/analysis/plot_v5.py',
        'reporting/analysis/plot_strengthening.py',
    ]
    for script in commands + plots:
        if script == plots[0]:
            for name in ['expansionrx_v5_identifiability_lead_v1', 'expansionrx_v6_strengthening_v1']:
                for source in (work / 'artifacts' / name).glob('*.csv'):
                    shutil.copyfile(source, work / 'reporting/tables' / source.name)
        log = work / 'logs' / (Path(script).stem + '.log')
        with log.open('wb') as stream:
            result = subprocess.run([sys.executable, '-B', '-X', 'utf8', script], cwd=work,
                                    stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f'{script} failed; see {log}')
        print(f'Completed {script}', flush=True)
    checks = []
    for entry in expected['files']:
        output = work / entry['path']
        if entry['kind'] == 'csv':
            passed = digest(output) == entry['sha256']
        else:
            actual = json.loads(output.read_text(encoding='utf-8'))
            for field in entry['excluded_paths']:
                parent = actual
                parts = field.split('.')
                for key in parts[:-1]:
                    parent = parent[key]
                parent.pop(parts[-1], None)
            passed = matches(actual, entry['content'], expected['json_comparison']['absolute_tolerance'])
        checks.append({'path': entry['path'], 'pass': passed})
    report = {'passed': sum(x['pass'] for x in checks), 'total': len(checks), 'checks': checks,
              'scope': 'Reanalysis of the supplied scores and selections; no model retraining or new predictions.'}
    (work / 'verification.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'passed': report['passed'], 'total': report['total'], 'output': str(work)}))
    if report['passed'] != report['total']:
        raise SystemExit('Reference-output verification failed; see verification.json')


if __name__ == '__main__':
    main()
