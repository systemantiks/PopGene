"""Packaging regression checks for staged dependency notice replacement."""
import importlib.util
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


if __name__ == '__main__':
    unittest.main()
