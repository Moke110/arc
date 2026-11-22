# HVG-DEG expression views

Read [README.md](README.md), ADRs 0062, 0082 and 0083, and the repository
execution environment before building and verifying views.

Use the bounded native-corpus and prediction-axis readers. Cap and monitor real
build and replay jobs, and run large workloads serially. Build in
`execution/checkpoints/` and verify there before publishing. Two-source
publication requires verified tenv and Tahoe corpus payloads and the fixed
5,000-gene axis. Preserve historical study payloads and sealed receipts.
