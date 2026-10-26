#!/usr/bin/env python3
"""Release safety checks; no network, signing credentials or app installation."""
import base64
import importlib.util
import json
import os
from pathlib import Path
import plistlib
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace("-", "_"), ROOT / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


package = load("package-macos")
release = load("prepare-release")
fetch = load("fetch-sparkle")


class ReleaseSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "release").mkdir()
        self.config = dict(build=12, feed_url="https://example.com/appcast.xml",
                           public_ed_key=base64.b64encode(bytes(range(32))).decode(),
                           keychain_account="test")
        self.write_config()

    def write_config(self):
        (self.root / "release/macos.json").write_text(json.dumps(self.config))

    def read_config(self, public=True):
        with patch.object(package, "ROOT", self.root), patch.dict(os.environ, {}, clear=True):
            return package.read_config(public)

    def test_release_requires_key_but_development_does_not(self):
        self.config["public_ed_key"] = ""
        self.write_config()
        self.assertEqual(self.read_config(False)["public_ed_key"], "")
        with self.assertRaisesRegex(ValueError, "before releasing"):
            self.read_config()

    def test_invalid_keys_and_insecure_feeds_fail_closed(self):
        for key in ["not base64", "AA=="]:
            self.config["public_ed_key"] = key
            self.write_config()
            with self.assertRaisesRegex(ValueError, "32-byte"):
                self.read_config()
        self.config["feed_url"] = "http://example.com/feed.xml"
        self.write_config()
        with self.assertRaisesRegex(ValueError, "HTTPS"):
            self.read_config()

    def test_plist_uses_monotonic_build_and_manual_signed_updates(self):
        plist = plistlib.loads(plistlib.dumps(package.bundle_plist("0.0.1", self.config)))
        self.assertEqual(plist["CFBundleShortVersionString"], "0.0.1")
        self.assertEqual(plist["CFBundleVersion"], "12")
        self.assertFalse(plist["SUEnableAutomaticChecks"])
        self.assertFalse(plist["SUAllowsAutomaticUpdates"])
        self.assertTrue(plist["SURequireSignedFeed"])
        self.assertTrue(plist["SUVerifyUpdateBeforeExtraction"])

    def archive(self, info):
        archive = self.root / "PopGene-0.0.1-macos-arm64.zip"
        with zipfile.ZipFile(archive, "w") as file:
            file.writestr("PopGene.app/Contents/Info.plist", plistlib.dumps(info))
        return archive

    def test_reject_stale_or_misconfigured_release_archive(self):
        info = package.bundle_plist("0.0.1", self.config)
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(release.validate_archive(self.archive(info), self.config), "0.0.1")
            for key, value in [("CFBundleVersion", "11"), ("SUPublicEDKey", ""),
                               ("SURequireSignedFeed", False), ("SUFeedURL", "https://wrong.example/feed"),
                               ("CFBundleIdentifier", "some.other.app")]:
                with self.subTest(key=key), self.assertRaises(ValueError):
                    release.validate_archive(self.archive({**info, key: value}), self.config)
        with patch.dict(os.environ, {"RELEASE_TAG": "v0.0.2"}):
            with self.assertRaisesRegex(ValueError, "tag"):
                release.validate_archive(self.archive(info), self.config)

    def test_build_numbers_cannot_go_backwards_or_be_reused(self):
        feed = self.root / "appcast.xml"
        feed.write_text('<rss xmlns:sparkle="' + release.SPARKLE + '"><channel><item>'
                        '<sparkle:version>12</sparkle:version></item></channel></rss>')
        for build in [11, 12]:
            with self.assertRaisesRegex(ValueError, "Increase"):
                release.validate_previous(feed, build)
        release.validate_previous(feed, 13)
        feed.write_text(feed.read_text().replace("12", "not-a-build"))
        with self.assertRaisesRegex(ValueError, "unsupported"):
            release.validate_previous(feed, 13)

    def test_framework_checksum_tampering_is_rejected(self):
        file = self.root / "framework.tar.xz"
        file.write_bytes(b"tampered")
        with self.assertRaisesRegex(SystemExit, "checksum mismatch"):
            fetch.verify(file, "0" * 64)


if __name__ == "__main__":
    unittest.main()
