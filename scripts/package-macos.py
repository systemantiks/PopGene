#!/usr/bin/env python3
"""Build a local app or signed, notarized macOS release with Sparkle."""
import argparse
import base64
import binascii
import importlib.util
import json
import os
from pathlib import Path
import plistlib
import re
import shutil
import subprocess
import tempfile
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]


def run(*args, **kwargs):
    return subprocess.run([str(arg) for arg in args], cwd=ROOT, check=True, **kwargs)


def read_config(release=False):
    config = json.loads((ROOT / "release/macos.json").read_text())
    config["public_ed_key"] = (os.environ.get("SPARKLE_PUBLIC_ED_KEY") or config["public_ed_key"]).strip()
    if type(config["build"]) is not int or config["build"] < 12:
        raise ValueError("Release build must be an increasing integer, starting at 12.")
    feed = urlparse(config["feed_url"])
    if feed.scheme != "https" or not feed.netloc or feed.username or feed.password:
        raise ValueError("The update feed must be an HTTPS URL without credentials.")
    key = config["public_ed_key"]
    if key:
        try:
            valid = len(base64.b64decode(key, validate=True)) == 32
        except (ValueError, binascii.Error):
            valid = False
        if not valid:
            raise ValueError("SPARKLE_PUBLIC_ED_KEY must be a base64 encoded 32-byte public key.")
    elif release:
        raise ValueError("Set public_ed_key in release/macos.json or SPARKLE_PUBLIC_ED_KEY before releasing.")
    return config


def bundle_plist(version, config):
    result = dict(
        CFBundleExecutable="PopGene", CFBundleIdentifier="org.popgene.desktop",
        CFBundleName="PopGene", CFBundleIconFile="PopGene.icns", CFBundlePackageType="APPL",
        CFBundleShortVersionString=version, CFBundleVersion=str(config["build"]),
        LSMinimumSystemVersion="12.0.0", NSHighResolutionCapable=True,
        SUFeedURL=config["feed_url"], SUEnableAutomaticChecks=False,
        SUAutomaticallyUpdate=False, SUAllowsAutomaticUpdates=False,
        SUSendProfileInfo=False, SUVerifyUpdateBeforeExtraction=True,
        SURequireSignedFeed=True,
    )
    if config["public_ed_key"]:
        result["SUPublicEDKey"] = config["public_ed_key"]
    return result


def sign_bundle(app, identity):
    framework = app / "Contents/Frameworks/Sparkle.framework"
    version = framework / "Versions/Current"
    # Sign inside out, retaining the Downloader helper's sandbox entitlement.
    targets = [version / "XPCServices/Downloader.xpc", version / "XPCServices/Installer.xpc",
               version / "Autoupdate", version / "Updater.app", framework, app]
    for target in targets:
        flags = ["--force", "--sign", identity]
        if identity != "-":
            flags += ["--options", "runtime", "--timestamp"]
        if target != app:
            flags += ["--preserve-metadata=identifier,entitlements"]
        run("codesign", *flags, target)
    run("codesign", "--verify", "--deep", "--strict", app)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", action="store_true", help="Require signing, notarization and cleared notices; create a ZIP")
    args = parser.parse_args()
    config = read_config(args.release)
    identity = os.environ.get("APPLE_SIGNING_IDENTITY", "-")
    if args.release and not identity.startswith("Developer ID Application:"):
        raise ValueError("Set APPLE_SIGNING_IDENTITY to a Developer ID Application identity.")
    profile = os.environ.get("APPLE_NOTARY_PROFILE")
    if args.release and not profile:
        raise ValueError("Set APPLE_NOTARY_PROFILE to your notarytool Keychain profile.")
    metadata = json.loads(run(ROOT / "scripts/cargo-local", "metadata", "--locked", "--offline",
                              "--no-deps", "--format-version=1", capture_output=True, text=True).stdout)
    version = next(p["version"] for p in metadata["packages"] if p["name"] == "popgene-desktop")
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ValueError("Stable macOS releases require a numeric major.minor.patch version.")
    if os.environ.get("RELEASE_TAG", "v" + version) != "v" + version:
        raise ValueError("Release tag must match the Cargo workspace version.")
    notices = ["python3", "scripts/collect-notices.py"]
    if args.release:
        notices.append("--strict")
    run(*notices)
    spec = importlib.util.spec_from_file_location("fetch_sparkle", ROOT / "scripts/fetch-sparkle.py")
    fetch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fetch)
    sparkle = fetch.fetch()
    run(ROOT / "scripts/cargo-local", "build", "--locked", "--release", "-p", "popgene-desktop")
    run("sh", "scripts/build-icons.sh")
    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="package-", dir=dist) as temp:
        app = Path(temp) / "PopGene.app"
        contents = app / "Contents"
        (contents / "MacOS").mkdir(parents=True)
        resources = contents / "Resources"
        resources.mkdir()
        shutil.copy2(Path(metadata["target_directory"]) / "release/popgene", contents / "MacOS/PopGene")
        shutil.copy2(ROOT / "LICENSE", resources / "LICENSE")
        shutil.copy2(ROOT / "assets/icons/PopGene.icns", resources / "PopGene.icns")
        shutil.copytree(dist / "third-party", resources / "third-party")
        shutil.copy2(sparkle / "LICENSE", resources / "third-party/Sparkle-LICENSE.txt")
        shutil.copy2(ROOT / "scripts/sparkle.json", resources / "third-party/Sparkle-version.json")
        (contents / "Frameworks").mkdir()
        run("ditto", sparkle / "Sparkle.framework", contents / "Frameworks/Sparkle.framework")
        with (contents / "Info.plist").open("wb") as file:
            plistlib.dump(bundle_plist(version, config), file)
        sign_bundle(app, identity)
        if args.release:
            upload = Path(temp) / "notarize.zip"
            run("ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", app, upload)
            notary_args = ["--keychain-profile", profile, "--wait"]
            if os.environ.get("APPLE_NOTARY_KEYCHAIN"):
                notary_args += ["--keychain", os.environ["APPLE_NOTARY_KEYCHAIN"]]
            run("xcrun", "notarytool", "submit", upload, *notary_args)
            run("xcrun", "stapler", "staple", app)
            run("xcrun", "stapler", "validate", app)
            run("spctl", "--assess", "--type", "execute", "--verbose=2", app)
            archive = dist / ("PopGene-" + version + "-macos-arm64.zip")
            architectures = run("lipo", "-archs", contents / "MacOS/PopGene", capture_output=True, text=True).stdout.strip()
            if architectures != "arm64":
                raise ValueError("This release pipeline currently supports macOS ARM64 only.")
            staged_archive = Path(temp) / archive.name
            # Sparkle must sign the final archive, after notarization and stapling.
            run("ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", app, staged_archive)
            staged_archive.replace(archive)
            print(archive)
        destination = dist / "PopGene.app"
        # Renaming preserves executable mappings held by an already-running app.
        old = Path(temp) / "Previous.app"
        if destination.exists():
            destination.rename(old)
        try:
            app.rename(destination)
        except OSError:
            if old.exists():
                old.rename(destination)
            raise
    print(destination)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))
