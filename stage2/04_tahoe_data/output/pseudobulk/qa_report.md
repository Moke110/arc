# Tahoe pseudobulk corpus — QA

- profiles: **67,018** (>=100 cells: **61,767**, 92.2%)
- median detected genes per profile: **25290** (gate >=12,000: PASS)
- LFC vs authors' DESeq2 DE over **3755** contrasts; genes with baseMean>800 (median 328/contrast): Pearson 0.812 (p10 0.669) - reported, not gated
- baseMean>500 (644/contrast): Pearson 0.779, Spearman 0.713 - gate >=0.60: PASS
- baseMean>300: Pearson 0.741; Spearman over all genes: 0.173 (reported, not gated - near-null genes dominate ranks)
- (plate, cell_line) pairs with a vehicle control: **700** checked, all present

Method: our log2((treated counts + 1) / (same-plate same-cell-line DMSO mean + 1)) vs DESeq2 `log2FoldChange`, per-contrast Pearson/Spearman restricted to the authors' well-measured genes (baseMean threshold). The full-gene Spearman is reported for context only: on genes with no real effect both estimates are noise, so rank agreement there is uninformative (measured ~0.2).
