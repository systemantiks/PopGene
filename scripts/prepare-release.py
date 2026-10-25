#!/usr/bin/env python3
"""Create signed GitHub release assets from the final notarized macOS ZIP."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import plistlib
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPARKLE = "http://www.andymatuschak.org/xml-namespaces/sparkle"


def validate_archive(archive, config):
    with zipfile.ZipFile(archive) as file:
        info = plistlib.loads(file.read("PopGene.app/Contents/Info.plist"))
    version = info["CFBundleShortVersionString"]
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ValueError("Only stable numeric versions can use the stable release feed.")
    if archive.name != "PopGene-" + version + "-macos-arm64.zip":
        raise ValueError("Archive name must include its exact version and architecture.")
    if info.get("CFBundleIdentifier") != "org.popgene.desktop":
        raise ValueError("Archive is not a PopGene application.")
    if info.get("CFBundleVersion") != str(config["build"]):
        raise ValueError("Archive build does not match release/macos.json; rebuild it.")
    if not info.get("SUPublicEDKey") or info["SUPublicEDKey"] != config["public_ed_key"]:
        raise ValueError("Archive's public key must match the configured release key.")
    if info.get("SUFeedURL") != config["feed_url"] or not info.get("SURequireSignedFeed"):
        raise ValueError("Archive must use the configured signed update feed.")
    if os.environ.get("RELEASE_TAG", "v" + version) != "v" + version:
        raise ValueError("Release tag does not match the archive version.")
    return version


def validate_previous(path, build):
    root = ET.parse(path).getroot()
    for item in root.findall("./channel/item"):
        version = item.findtext("{" + SPARKLE + "}version")
        if version is None:
            enclosure = item.find("enclosure")
            version = enclosure.get("{" + SPARKLE + "}version") if enclosure is not None else None
        if version is None or not version.isdigit():
            raise ValueError("Previous feed contains an unsupported build number.")
        if int(version) >= build:
            raise ValueError("Increase release/macos.json build above every published build.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--notes", required=True, type=Path, help="Markdown release notes")
    parser.add_argument("--previous-appcast", type=Path, help="Required after the first published release")
    parser.add_argument("--output", type=Path, default=ROOT / "dist/release")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Output already exists; choose a fresh --output directory.")
    spec = importlib.util.spec_from_file_location("package_macos", ROOT / "scripts/package-macos.py")
    package = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(package)
    config = package.read_config(release=True)
    version = validate_archive(args.archive, config)
    notes = args.notes.read_text()
    if not notes.strip():
        raise ValueError("Release notes must not be empty.")
    spec = importlib.util.spec_from_file_location("fetch_sparkle", ROOT / "scripts/fetch-sparkle.py")
    fetch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fetch)
    sparkle = fetch.fetch()
    secret = os.environ.get("SPARKLE_PRIVATE_KEY")
    key_args = ["--ed-key-file", "-"] if secret else ["--account", config["keychain_account"]]

    def signing_tool(name, *arguments):
        # Never put the private key in argv, a log, or a checked-in file.
        subprocess.run([str(sparkle / "bin" / name), *key_args, *map(str, arguments)],
                       input=secret + "\n" if secret else None, text=True, check=True)

    if args.previous_appcast:
        signing_tool("sign_update", "--verify", args.previous_appcast)
        validate_previous(args.previous_appcast, config["build"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="release-", dir=args.output.parent) as work:
        staged = Path(work) / "assets"
        staged.mkdir()
        shutil.copy2(args.archive, staged / args.archive.name)
        (staged / (args.archive.stem + ".md")).write_text(notes)
        if args.previous_appcast:
            shutil.copy2(args.previous_appcast, staged / "appcast.xml")
        signing_tool("generate_appcast", "--maximum-deltas", "0", "--maximum-versions", "0",
                     "--embed-release-notes", "--download-url-prefix",
                     "https://github.com/systemantiks/PopGene/releases/download/v" + version + "/",
                     "--link", "https://github.com/systemantiks/PopGene/releases",
                     "-o", staged / "appcast.xml", staged)
        signing_tool("sign_update", "--verify", staged / "appcast.xml")
        # The upstream generator verifies the signing key against the archived
        # app's SUPublicEDKey. Confirm it actually emitted our build and signature.
        root = ET.parse(staged / "appcast.xml").getroot()
        items = [item for item in root.findall("./channel/item")
                 if item.findtext("{" + SPARKLE + "}version") == str(config["build"])]
        enclosure = items[0].find("enclosure") if len(items) == 1 else None
        if enclosure is None or not enclosure.get("{" + SPARKLE + "}edSignature"):
            raise ValueError("Sparkle did not produce a signed update for this build.")
        # Keep release notes for the GitHub release body as well as the feed.
        (staged / "release-notes.md").write_text(notes)
        (staged / (args.archive.stem + ".md")).unlink()
        if args.output.exists():
            raise ValueError("Output already exists; choose a fresh --output directory.")
        staged.rename(args.output)
    print(args.output)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))
