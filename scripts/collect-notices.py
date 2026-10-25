"""Collect target-specific dependency notices without retaining stale packages.

All downloads are reviewed and committed under vendor-notices before packaging;
this collector itself is offline. --strict also rejects recorded provenance gaps.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def is_notice(path):
    """Recognize notice text, not source files such as objc2's copying.rs."""
    if not path.is_file() or path.suffix.lower() in {
        '.rs', '.py', '.c', '.cpp', '.h', '.js', '.ts', '.astro', '.html', '.svg'
    }:
        return False
    name = path.name.upper()
    return (name.startswith(('LICENSE', 'LICENCE', 'COPYING', 'NOTICE', 'COPYRIGHT'))
            or name in ('OFL.TXT', 'UFL.TXT', 'UNLICENSE')
            or ('LICENSE' in name and path.suffix.lower() in ('.txt', '.md', '.rst')))


def replace_directory(staged, destination):
    """Keep the previous complete inventory until the new tree is ready."""
    if destination.is_symlink():
        raise RuntimeError(f'Refusing to replace a symlink: {destination}')
    previous = staged.parent / 'previous'
    had_previous = destination.exists()
    if had_previous:
        destination.rename(previous)
    try:
        staged.rename(destination)
    except BaseException:
        if had_previous:
            previous.rename(destination)
        raise
    if had_previous:
        shutil.rmtree(previous)


def collect(metadata, target, destination, strict=False):
    resolved = {node['id'] for node in metadata['resolve']['nodes']}
    supplements = json.loads((ROOT / 'vendor-notices/manifest.json').read_text())
    gaps_path = ROOT / 'vendor-notices/provenance-gaps.json'
    gaps = json.loads(gaps_path.read_text()) if gaps_path.exists() else {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    inventory = []
    missing = []
    with tempfile.TemporaryDirectory(prefix='.third-party-stage-', dir=destination.parent) as temporary:
        staged = Path(temporary) / 'notices'
        staged.mkdir()
        for package in sorted(metadata['packages'], key=lambda p: (p['name'], p['version'])):
            if package['id'] not in resolved or not package['source']:
                continue
            source = Path(package['manifest_path']).parent
            key = f"{package['name']}-{package['version']}"
            out = staged / key
            notices = sorted(p for p in source.rglob('*') if is_notice(p))
            license_file = package.get('license_file')
            if license_file:
                explicit = source / license_file
                if explicit.is_file() and explicit not in notices:
                    notices.append(explicit)
            for file in notices:
                # Registry manifests must not make the collector copy outside a crate.
                relative = file.resolve().relative_to(source.resolve())
                output = out / relative
                output.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(file, output)
            supplement = supplements.get(key)
            supplement_files = []
            if supplement:
                supplement_source = ROOT / 'vendor-notices' / supplement
                supplement_source.resolve().relative_to((ROOT / 'vendor-notices').resolve())
                if not supplement_source.is_dir():
                    raise RuntimeError(f'Missing supplement directory for {key}: {supplement}')
                supplement_files = sorted(str(p.relative_to(supplement_source)) for p in supplement_source.rglob('*') if p.is_file())
                shutil.copytree(supplement_source, out / 'upstream-supplement')
            if not notices and not supplement_files:
                missing.append(key)
            provenance = gaps.get(key)
            entry = {
                'name': package['name'], 'version': package['version'],
                'license': package['license'], 'authors': package.get('authors', []),
                'repository': package['repository'],
                'notice_files': len(notices),
                'notice_paths': [str(p.resolve().relative_to(source.resolve())) for p in notices],
                'supplement': supplement, 'supplement_files': supplement_files,
                'provenance_gap': provenance,
            }
            inventory.append(entry)
        current_gaps = {f"{p['name']}-{p['version']}": p['provenance_gap'] for p in inventory if p['provenance_gap']}
        if missing:
            raise RuntimeError('No notices or reviewed supplements for: ' + ', '.join(missing))
        if strict and current_gaps:
            raise RuntimeError('Strict notice review requires resolving: ' + ', '.join(current_gaps))
        result = {
            'target': target,
            'cargo_lock_sha256': hashlib.sha256((ROOT / 'Cargo.lock').read_bytes()).hexdigest(),
            'scope': 'Resolved registry dependencies for this target, including build dependencies; not a linked-code SBOM.',
            'provenance_gaps': current_gaps,
            'packages': inventory,
        }
        (staged / 'inventory.json').write_text(json.dumps(result, indent=2) + '\n')
        shutil.copytree(ROOT / 'vendor-notices/biopython-reference', staged / 'biopython-format-reference')
        shutil.copy2(ROOT / 'THIRD_PARTY.md', staged / 'README.md')
        replace_directory(staged, destination)
    print(f'Collected notices and reviewed supplements for {len(inventory)} resolved packages ({target}).')
    print('Packages without notices or supplements: none.')
    if current_gaps:
        print('Recorded upstream provenance gaps: ' + ', '.join(current_gaps))
        print('These are explicit in inventory.json; --strict rejects them for distribution review.')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', help='Rust target triple; defaults to the toolchain host')
    parser.add_argument('--strict', action='store_true', help='Fail on recorded upstream provenance gaps too')
    args = parser.parse_args()
    cargo = str(ROOT / 'scripts/cargo-local')
    target = args.target
    if target is None:
        version = subprocess.check_output([cargo, '-vV'], cwd=ROOT, text=True)
        target = next(line.split(': ', 1)[1] for line in version.splitlines() if line.startswith('host: '))
    metadata = json.loads(subprocess.check_output([
        cargo, 'metadata', '--locked', '--offline', '--format-version', '1', '--filter-platform', target
    ], cwd=ROOT, text=True))
    collect(metadata, target, ROOT / 'dist/third-party', args.strict)


if __name__ == '__main__':
    main()
