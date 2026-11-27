# Select UniMol for subsequent CarPS-VAE studies

**Status** Accepted, 2026-11-27.

Close the [Step 08 chemical representation comparison](../../stage2_CarPS-VAE/08_cpsv_chem-emb/README.md#formal-results) and use its frozen UniMol representation in subsequent CarPS-VAE studies, starting with the Step 09 gene-axis comparison.

Across nine matched cold-drug fold and training-seed coordinates, UniMol had the highest mean held-out Δ PCC (0.441262) and the lowest mean held-out Δ RMSE (0.238599) among the seven representations. Relative to KPGT, the paired changes were +0.012432 ± 0.013493 in Δ PCC and −0.002153 ± 0.002016 in Δ RMSE (mean ± sample SD). The mean Δ PCC gain was smaller than its sample SD. This successor configuration adopts UniMol by project decision under the separate adoption route in [ADR 0099](0099-carps-vae-alpha-and-chemical-representations.md).
