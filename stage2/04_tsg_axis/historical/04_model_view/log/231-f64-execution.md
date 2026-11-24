# Approved float64 formal execution

2026-11-17. User approved “批准 float64 变体并继续 #226”; source https://github.com/Moke110/arc/issues/226#issuecomment-6936759990 and https://app.notion.com/p/41bab6347e96818fafa7e11ffe71dec7.

Reviewed runtime9975fc69ac016cbc7ac7cc323834798a488daa09 (same tree as reviewed d9276fd). Scientific and runtime code identities are frozen in 231-f64-frozen-identity.json. Compatibility audit is prnet-f64-compatibility-audit.json. Candidate: /tmp/arc-stage2-tox-design/stage2/02_models_tenv-5k/execution/checkpoints/matrix-f64-v1. Original runtime/candidate retained unchanged. Inputs /tmp/arc-tenv5k-inputs and means /tmp/arc-tenv5k-means; source trees /home/chang/arc/pert_model/{PRnet,TranSiGen}/source.

Commands are recorded in /tmp/arc-tenv5k-resolve/run-formal-f64.sh. All heavy phases use execution/capped.sh, /usr/bin/time -v, CUDA native escalation and two threads. Logs/time evidence /tmp/arc-tenv5k-resolve/formal-f64-{import,profile,matrix}.{stdout,stderr,time}. Coordinator monitors shared GPU/WSL resources every10s. Import15 with original-runtime fresh replay, profile one formal identity then resume75 new runs unchanged. No publishing authorization delegated to this executor.
