# Approved float64 formal execution

2026-11-17. User approved “批准 float64 变体并继续 #226”; source https://github.com/Moke110/arc/issues/226#issuecomment-6936759990 and https://app.notion.com/p/41bab6347e96818fafa7e11ffe71dec7.

Reviewed runtimef5b8909ca8c9c25daa0912b582a62ecaaa50d9bc (final runtime source enumeration fix c263c2f merged; three PRnet method dependencies unchanged from reviewed d9276fd). Scientific and runtime code identities are frozen in 231-f64-v2-frozen-identity.json. Compatibility audit is prnet-f64-v2-compatibility-audit.json. Candidate: /tmp/arc-stage2-tox-design/stage2/02_models_tenv-5k/execution/checkpoints/matrix-f64-v2. Original runtime/candidate retained unchanged. Inputs /tmp/arc-tenv5k-inputs and means /tmp/arc-tenv5k-means; source trees /home/chang/arc/pert_model/{PRnet,TranSiGen}/source.

Commands are recorded in /tmp/arc-tenv5k-resolve/run-formal-f64-v2.sh. All heavy phases use execution/capped.sh, /usr/bin/time -v, CUDA native escalation and two threads. Logs/time evidence /tmp/arc-tenv5k-resolve/formal-f64-v2-{import,profile,matrix}.{stdout,stderr,time}. Coordinator monitors shared GPU/WSL resources every10s. Import15 with original-runtime fresh replay, profile one formal identity then resume75 new runs unchanged. No publishing authorization delegated to this executor.

Preserved matrix-f64-v1 unchanged after its engineering identity blocker. Runtime enumeration fix reviews: /tmp/arc-tenv5k-resolve/review-runtime-hash-{standards,spec}.md; merge receipt runtime-hash-merge.md.
