# macOS releases and updates

PopGene uses Sparkle 2.10.0 for **PopGene → Check for Updates…**. Checks are manual;
opening the app does not contact an update server. Sparkle displays release notes,
downloads and verifies the archive, installs it and relaunches the application.
PopGene postpones the restart for its existing save/discard/cancel flow. Cancelling
keeps the app open; choose Check for Updates again when ready to restart. In-progress
file operations or open dialogs defer the restart too.

The stable feed is an asset of the latest published GitHub release:

`https://github.com/systemantiks/PopGene/releases/latest/download/appcast.xml`

Every release must include both the versioned ZIP and `appcast.xml`. Keep the
repository and release assets publicly downloadable. No GitHub token is embedded
in the app. Draft releases and GitHub prereleases are not the stable update
channel. There is no need to enable GitHub Pages. This pipeline currently builds
for **Apple silicon only**; Windows/Linux self-updating is not implemented.

## Version and trust configuration

- `Cargo.toml` is the source of the displayed version, currently `0.0.1`.
- `release/macos.json` stores the build counter (starting at 12, after the old
  development build 11), feed URL, public verification key and Keychain account.
- Increment the build counter for every published binary, even if its displayed
  version stays the same. Never reuse a published tag, versioned ZIP or build.
- The public key is safe to commit. The corresponding private signing key is in
  the macOS login Keychain under account `org.popgene.desktop`, service
  `https://sparkle-project.org`, created during initial updater setup.
- `scripts/sparkle.json` pins the upstream distribution and its SHA-256. Both
  archives and the feed are signed; downloads are verified before extraction.

Back up the private key securely before distributing the first release. Losing it
can prevent existing installations from accepting future updates. Do not generate
a replacement key for each build. Key rotation requires a separate migration.

## Local builds

```sh
scripts/package-macos.sh
open dist/PopGene.app
```

This creates an ad-hoc signed app containing Sparkle and its notices. `cargo run`
also builds without Sparkle installed, but its unbundled executable cannot update
itself. Before the first release is published, the configured feed will return a
network error; it is not an indication that the development build is up to date.

## One-time distribution setup

1. Install a **Developer ID Application** signing certificate and private key
   from your Apple Developer account. A development certificate is insufficient
   for distribution. Store notarization credentials with `xcrun notarytool
   store-credentials popgene-release` and follow its prompts.
2. Resolve the existing notice provenance gaps listed in `THIRD_PARTY.md`.
   Public packaging deliberately runs `collect-notices.py --strict` and stops
   while any remain. Local development packaging continues to record the gaps.
3. On a new signing Mac, restore the existing Sparkle key. The pinned
   `bin/generate_keys --account org.popgene.desktop -p` prints only its public key.
   `-x /secure/path/private-key` exports a backup; `-f /secure/path/private-key`
   imports it. Protect that file and never add it to the repository.

The release scripts use the Keychain by default. In CI, `SPARKLE_PRIVATE_KEY` is
passed to Sparkle over stdin, never as a command-line argument or checked-in file.
`SPARKLE_PUBLIC_ED_KEY` can override the checked-in public key, but must match the
private key and the already-shipped app. Leave it unset for the default key.

## Build a release locally

Update the workspace version, run Cargo once to update `Cargo.lock`, increment
`release/macos.json`'s build, write `release/notes/<version>.md` and update the
versioned download links in `README.md`. Then run:

```sh
export APPLE_SIGNING_IDENTITY='Developer ID Application: Your Name (TEAMID)'
export APPLE_NOTARY_PROFILE=popgene-release
scripts/package-macos.sh --release

# First release:
python3 scripts/prepare-release.py \
  --archive dist/PopGene-0.0.1-macos-arm64.zip \
  --notes release/notes/0.0.1.md
```

The packaging step signs helpers inside out, notarizes, staples, verifies and
archives the final app. The feed generator signs those exact bytes and verifies
the generated feed. It refuses stale build numbers, missing keys, mismatched tags
or an existing output directory. Move an old `dist/release` aside before a new run.

For subsequent releases, download the previous stable release's `appcast.xml`
using `gh release download <previous-tag> --pattern appcast.xml --dir dist/previous`
and pass `--previous-appcast dist/previous/appcast.xml` to `prepare-release.py`.
It verifies the old feed before extending it, preserves historical download URLs
and requires the new build to exceed all published builds. This supports older
macOS versions when a future release raises its minimum OS requirement.

## GitHub Actions

`.github/workflows/release-macos.yml` runs for pushed `v*` tags or an explicitly
selected existing tag via **Actions → Draft macOS release → Run workflow**.
It validates, packages and notarizes, downloads the previous feed, signs assets,
and creates a **draft** GitHub release. It does not publish automatically.

Configure the GitHub environment named `release` with:

| Type | Name | Value |
| --- | --- | --- |
| Variable | `APPLE_SIGNING_IDENTITY` | Full Developer ID Application identity |
| Variable, optional | `SPARKLE_PUBLIC_ED_KEY` | Override only if using the corresponding signing key |
| Secret | `APPLE_CERTIFICATE_P12_BASE64` | Base64 of the exported signing certificate plus private key |
| Secret | `APPLE_CERTIFICATE_PASSWORD` | Password protecting that P12 |
| Secret | `APPLE_ID` | Apple account used for notarization |
| Secret | `APPLE_TEAM_ID` | Apple Developer team ID |
| Secret | `APPLE_APP_PASSWORD` | App-specific password for notarization |
| Secret | `SPARKLE_PRIVATE_KEY` | Contents of the Sparkle key backup export |

The workflow creates a temporary signing Keychain and deletes it even on failure.
For repeat builds after failure, inspect any existing draft before deleting or
replacing it; the workflow intentionally does not overwrite release assets.

After reviewing the draft's ZIP, feed and notes, publish it as the **latest stable
release**. The first installation is a manual download; subsequent versions can
update from the app. Never mark an incomplete release as latest. When creating a
release through the GitHub UI, attach the prepared ZIP and feed before publishing.

## Verification

```sh
python3 scripts/test-release.py
scripts/cargo-local test --locked --workspace
scripts/cargo-local clippy --locked --workspace --all-targets -- -D warnings
scripts/cargo-local fmt --all --check

# Native termination regression, in a graphical macOS session:
clang -fobjc-arc -Wall -Wextra -Werror -framework AppKit \
  crates/popgene-desktop/src/updater_macos.m \
  validation/macos-updater-termination.m -o /tmp/popgene-termination-test
/tmp/popgene-termination-test
```

Before public distribution, exercise an update between two signed/notarized test
builds: no update, new update, offline check, invalid signature, Save/Discard/Cancel,
save failure, busy file I/O, install, relaunch and version display. Use a separate
test feed and test bundle identity so test versions never reach the public feed.
The native regression specifically checks cancellation and retry of a postponed
Sparkle restart without entering AppKit's modal quit loop, which blocks GPUI.

Local validation on 23 October 2026 passed 158 Rust tests, six release safety
tests, clippy, formatting, the native termination regression and a signed release
asset generation smoke test. An isolated GPUI app copy downloaded a signed update
from localhost, showed its unsaved-changes prompt, cancelled and retried, installed
build 13 over build 12, automatically relaunched, and reported itself up to date.
A modified feed was rejected with a signature error. These trial bundles used
the same application executable with different bundle versions and ad-hoc code
signatures; Developer ID signing/notarization and GitHub Actions still require
the distribution credentials and notice review above.

Upstream references: [Sparkle setup](https://sparkle-project.org/documentation/),
[publishing updates](https://sparkle-project.org/documentation/publishing/),
[GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository).
