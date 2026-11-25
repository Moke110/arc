# Condition the tenv axis comparison on matched plate controls

**Status:** accepted (2026-11-25). See the [formal rerun receipt](../../stage2/04_tsg_axis/log/319-formal.json).

The Step 04 successor measures response prediction given a measured, matched same-plate DMSO profile. For each condition, use the Step 03 plate control on its prediction axis, `C_b`, as the control expression input in training, validation and test. The model predicts treated expression. Compute predicted and measured Delta from the same control, `predicted Delta = predicted treatment − C_b` and `measured Delta = measured treatment − C_b`. The original [Spec 5](https://github.com/Moke110/arc/issues/310) protocol used the axis-wide reference for held-out inference and predicted Delta, which measured performance under a different control-availability assumption.

Use the author-published TranSiGen MSE + KLD objective with `beta = 0.1` for training and validation. Choose the checkpoint with the lowest per-condition mean validation objective and stop after 30 consecutive epochs without improvement. Replay validation and held-out predictions from the selected checkpoint. Keep test measurements for scoring.

Fit all six axes on all three frozen cold-drug folds from scratch. Publish the successor metrics and figures after acceptance and preserve the original formal release and its pinned replay.
