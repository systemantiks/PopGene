"""Packaging regression checks for staged dependency notice replacement."""
import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('collect_notices', Path(__file__).with_name('collect-notices.py'))
collector = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(collector)


class CollectorTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.original_root = collector.ROOT
        collector.ROOT = self.root
        self.addCleanup(setattr, collector, 'ROOT', self.original_root)
        vendor = self.root / 'vendor-notices'
        (vendor / 'biopython-reference').mkdir(parents=True)
        (vendor / 'biopython-reference/LICENSE.rst').write_text('reference licence')
        for relative in ['crates/popgene-core/vendor/primer3', 'vendor-notices/sgffp-reference',
                         'vendor-notices/autosnapgene-reference',
                         'vendor-notices/staden-io-lib', 'vendor-notices/figure-fonts']:
            component = self.root / relative
            component.mkdir(parents=True)
            (component / 'LICENSE.txt').write_text('local component licence')
        (vendor / 'manifest.json').write_text('{}')
        (self.root / 'THIRD_PARTY.md').write_text('Reviewed notices')
        (self.root / 'Cargo.lock').write_text('lock')
        self.source = self.root / 'registry/example-1.0.0'
        self.source.mkdir(parents=True)
        (self.source / 'Cargo.toml').write_text('[package]')
        self.destination = self.root / 'dist/third-party'
        self.destination.mkdir(parents=True)
        (self.destination / 'stale-egui-notice').write_text('previous')
        self.metadata = {
            'resolve': {'nodes': [{'id': 'example'}]},
            'packages': [{'id': 'example', 'name': 'example', 'version': '1.0.0',
                          'manifest_path': str(self.source / 'Cargo.toml'),
                          'source': 'registry', 'license': 'MIT', 'authors': ['Author'],
                          'repository': 'https://example.invalid/project'}],
        }

    def test_replacement_removes_stale_packages_and_preserves_actual_notices(self):
        (self.source / 'LICENSE-MIT').write_text('licence text')
        (self.source / 'copying.rs').write_text('not a licence')
        result = collector.collect(self.metadata, 'test-target', self.destination)
        self.assertFalse((self.destination / 'stale-egui-notice').exists())
        self.assertEqual((self.destination / 'example-1.0.0/LICENSE-MIT').read_text(), 'licence text')
        self.assertFalse((self.destination / 'example-1.0.0/copying.rs').exists())
        self.assertEqual(result['packages'][0]['notice_files'], 1)

    def test_missing_notice_does_not_replace_existing_inventory(self):
        with self.assertRaisesRegex(RuntimeError, 'No notices'):
            collector.collect(self.metadata, 'test-target', self.destination)
        self.assertEqual((self.destination / 'stale-egui-notice').read_text(), 'previous')

    def test_strict_provenance_check_leaves_previous_inventory_intact(self):
        (self.source / 'LICENSE').write_text('licence text')
        (self.root / 'vendor-notices/provenance-gaps.json').write_text(json.dumps({
            'example-1.0.0': {'kind': 'later_upstream_notice'}
        }))
        with self.assertRaisesRegex(RuntimeError, 'Strict notice review'):
            collector.collect(self.metadata, 'test-target', self.destination, strict=True)
        self.assertTrue((self.destination / 'stale-egui-notice').exists())

    def test_strict_accepts_reviewed_provenance_gap_and_keeps_it_in_inventory(self):
        (self.source / 'LICENSE').write_text('licence text')
        (self.root / 'vendor-notices/provenance-gaps.json').write_text(json.dumps({
            'example-1.0.0': {'kind': 'later_upstream_notice',
                              'review': {'status': 'accepted', 'date': '2026-10-26', 'basis': 'reviewed'}}
        }))
        result = collector.collect(self.metadata, 'test-target', self.destination, strict=True)
        self.assertIn('example-1.0.0', result['provenance_gaps'])
        self.assertFalse((self.destination / 'stale-egui-notice').exists())

    def test_failed_install_restores_previous_directory(self):
        stage = self.root / 'stage'
        stage.mkdir()
        staged = stage / 'notices'
        staged.mkdir()
        original = Path.rename

        def fail_new_tree(path, target):
            if path == staged:
                raise OSError('simulated rename failure')
            return original(path, target)

        with patch.object(Path, 'rename', fail_new_tree):
            with self.assertRaisesRegex(OSError, 'simulated'):
                collector.replace_directory(staged, self.destination)
        self.assertEqual((self.destination / 'stale-egui-notice').read_text(), 'previous')

    def test_embedded_font_licenses_and_checked_assets_are_included(self):
        with patch.object(collector, 'ROOT', self.original_root):
            stage = self.root / 'real-notices'
            stage.mkdir()
            components = collector.collect_local_components(stage)
        fonts = next(item for item in components if item['name'] == 'figure-fonts')
        self.assertEqual(fonts['source'], 'vendor-notices/figure-fonts')
        provenance = json.loads((stage / 'figure-fonts/SOURCE.json').read_text())
        self.assertEqual(len(provenance['files']), 2)
        for entry in provenance['files']:
            data = (self.original_root / entry['asset']).read_bytes()
            self.assertEqual(hashlib.sha256(data).hexdigest(), entry['sha256'])
            license_data = (stage / 'figure-fonts' / entry['license_path']).read_bytes()
            self.assertEqual(hashlib.sha256(license_data).hexdigest(), entry['license_sha256'])
            self.assertIn(b'SIL OPEN FONT LICENSE Version 1.1', license_data)
            self.assertEqual(fonts['sha256'][entry['license_path']], entry['license_sha256'])


if __name__ == '__main__':
    unittest.main()
