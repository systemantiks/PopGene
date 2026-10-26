#!/bin/sh
# Derive the macOS icon family from the checked-in master artwork.
set -eu
ROOT=$(CDPATH='' cd -- "$(dirname -- "$0")/.." && pwd)
ICONS="$ROOT/assets/icons"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/popgene-icons.XXXXXX")
trap 'rm -rf "$WORK"' EXIT HUP INT TERM
mkdir -p "$WORK/PopGene.iconset"
for SIZE in 16 32 128 256 512; do
    sips -z "$SIZE" "$SIZE" "$ICONS/popgene-source.png" \
        --out "$WORK/PopGene.iconset/icon_${SIZE}x${SIZE}.png" >/dev/null
    DOUBLE=$((SIZE * 2))
    sips -z "$DOUBLE" "$DOUBLE" "$ICONS/popgene-source.png" \
        --out "$WORK/PopGene.iconset/icon_${SIZE}x${SIZE}@2x.png" >/dev/null
done
iconutil -c icns "$WORK/PopGene.iconset" -o "$ICONS/PopGene.icns"
printf '%s\n' "$ICONS/PopGene.icns"
