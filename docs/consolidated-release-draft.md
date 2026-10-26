# Consolidated local build — 0.0.1

The consolidated changes and current validation are recorded in [local build notes](../release/notes/0.0.1-consolidated.md). These replace the earlier chronological release draft; they are not a public release announcement or a claim of complete parity.

The source0017 baseline passed 443 workspace tests, formatting, strict lint and debug linking on macOS ARM64 and Linux ARM64. Windows GNU compile/link checks passed without runtime verification. A subsequent narrow circular-RNA import correction passed 172 formats tests, scoped strict lint, formatting and a byte-preserving reference-file probe; the complete platform matrix was not repeated after that patch.

Focused macOS GUI checks now cover crash recovery, GenBank/EMBL physical-fragment save/close/reopen/rejoin, the exact native rejoined product in the reference application, explicit Grimsley selection and result invalidation, and stable uncolored-primer display after native Save/close/reopen. The optimized app also passed Open cancellation/retry, multi-selection and subsequent opens after restarting the earlier stalled temporary app. That historical chooser stall has no established cause or claimed code fix.

Final post-patch artifact identity, installation and last smoke-check results are recorded by the release owner in the local build notes. Format, scientific-model and platform limits remain explicit in [format support](file-formats.md), [protein models](protein-models.md), [validation history](development.md) and the [implementation checklist](implementation-checklist.md).
