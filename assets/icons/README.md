# PopGene app icon

The current app icon is an orange cat holding a teal-and-pink DNA strand.

- `popgene-source.png`: master artwork (1254 × 1254, RGBA).
- `PopGene.icns`: macOS icon family, 16–1024 pixels including Retina representations.

Run `sh scripts/build-icons.sh` on macOS to regenerate the icon family using Apple's `sips` and `iconutil`. Temporary size variants are removed after packaging. `scripts/package-macos.sh` includes the icon in the app bundle and declares it through `CFBundleIconFile` before signing.
