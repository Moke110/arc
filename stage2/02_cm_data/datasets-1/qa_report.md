# datasets-1 — QA report

Corpus: `corpus/v1`

Gene axis: **18,164** genes

Compounds in corpus: **525** drug keys

## Dataset ids

| dataset_id | accession |
|---|---|
| `dtoxs_primary` | — |
| `dtoxs_ipsc` | — |
| `screen_ctox` | GSE244740 |
| `temp_env` | GSE262419 |

## Per-view counts

| view | perturbed | controls | drug keys | obs |
|---|---|---|---|---|
| `dtoxs_cm.h5ad` | 957 | 118 | 46 | 1,075 |
| `screen_ctox.h5ad` | 2,520 | 552 | 42 | 3,072 |
| `temp_env_atc-moa.h5ad` | 1,047 | 272 | 113 | 1,319 |
| `temp_env_atc-moa-mesh.h5ad` | 2,310 | 272 | 252 | 2,582 |
| `temp_env_all-smiles.h5ad` | 4,206 | 272 | 460 | 4,478 |

## Gates

| gate | result | detail |
|---|---|---|
| var identical across all views | PASS | 18,164 genes |
| tier nesting temp_env_atc-moa.h5ad ⊆ temp_env_atc-moa-mesh.h5ad | PASS | temp_env_atc-moa.h5ad rows absent from temp_env_atc-moa-mesh.h5ad: [] |
| tier nesting temp_env_atc-moa-mesh.h5ad ⊆ temp_env_all-smiles.h5ad | PASS | temp_env_atc-moa-mesh.h5ad rows absent from temp_env_all-smiles.h5ad: [] |
| temp_env controls identical across views | PASS | 272 control rows |
| tier counts dtoxs | PASS | {'atc-moa': 46, 'atc-moa-mesh': 46, 'all-smiles': 47} |
| tier counts screen_ctox | PASS | {'atc-moa': 42, 'atc-moa-mesh': 42, 'all-smiles': 42} |
| tier counts temp_env | PASS | {'atc-moa': 113, 'atc-moa-mesh': 252, 'all-smiles': 460} |
| source values dtoxs_primary | PASS | 2 sampled cells match the source file |
| source values dtoxs_ipsc | PASS | 2 sampled cells match the source file |
| source values screen_ctox | PASS | 2 sampled cells match the source file |
| source values temp_env | PASS | 2 sampled cells match the source file |
| stage1 dtoxs_primary: well set | PASS | 662 Stage 1 wells all present; 180 extra source rows (controls / HG-excluded) |
| stage1 dtoxs_primary: doses | PASS | 662 shared wells |
| stage1 dtoxs_primary: counts rank | PASS | 962 shared genes, 5 sampled wells, min Spearman 1.000000 |
| stage1 dtoxs_ipsc: well set | PASS | 296 Stage 1 wells all present; 48 extra source rows (controls / HG-excluded) |
| stage1 dtoxs_ipsc: doses | PASS | 296 shared wells |
| stage1 dtoxs_ipsc: counts rank | PASS | 962 shared genes, 5 sampled wells, min Spearman 1.000000 |
| stage1 screen_ctox: well set | PASS | 2,424 Stage 1 wells all present; 648 extra source rows (controls / HG-excluded) |
| stage1 screen_ctox: doses | PASS | 2,424 shared wells; 48 wells corrected from the Stage 1 underscore-name fallback (['Amphotericin B', 'Mitomycin C']) |
| stage1 screen_ctox: counts rank | PASS | 962 shared genes, 5 sampled wells, min Spearman 1.000000 |
| stage1 temp_env: well set | PASS | 936 Stage 1 wells all present; 3,622 extra source rows (controls / HG-excluded) |
| stage1 temp_env: doses | PASS | 936 shared wells |
| stage1 temp_env: counts rank | PASS | 638 shared genes, 5 sampled wells, min Spearman 1.000000; 290 Stage 1 entrez-mapped landmarks excluded (probe-id suffix bug, recovered by symbol); 35 probe-suffix collision genes excluded (Stage 1 hijacks the probe by suffix) |

## Gene-side losses

Per source: identifiers seen, how they resolved, and the unique genes
they contributed. `unmapped`/`ambiguous` identifiers are dropped, never
guessed.

| source | direct | alias | previous | entrez | ensembl | entrez-verified | unmapped | ambiguous | unique genes | in intersection |
|---|---|---|---|---|---|---|---|---|---|---|
| dtoxs_primary | 23,410 | 70 | 877 | 2,047 | 0 | 25,445 | 12,073 | 1 | 26,328 | 18,164 |
| dtoxs_ipsc | 23,410 | 70 | 877 | 2,047 | 0 | 25,445 | 12,073 | 1 | 26,328 | 18,164 |
| screen_ctox | 0 | 0 | 0 | 0 | 31,880 | 0 | 12,274 | 0 | 31,880 | 18,164 |
| temp_env | 21,733 | 3 | 705 | 17 | 0 | 2 | 79 | 0 | 19,613 | 18,164 |

## Excluded compounds (not in the corpus)

| dataset | compound | entity | wells | reason |
|---|---|---|---|---|
| dtoxs_ipsc | Bevacizumab | antibody | 8 | non-small-molecule perturbagen (antibody) |
| dtoxs_ipsc | Cetuximab | antibody | 8 | non-small-molecule perturbagen (antibody) |
| dtoxs_ipsc | Rituximab | antibody | 8 | non-small-molecule perturbagen (antibody) |
| dtoxs_ipsc | Trastuzumab | antibody | 8 | non-small-molecule perturbagen (antibody) |
| dtoxs_primary | Endothelin-1 | peptide | 14 | non-small-molecule perturbagen (peptide) |
| dtoxs_primary | IGF-1 | protein | 16 | non-small-molecule perturbagen (protein) |
| dtoxs_primary | TNF-alpha | protein | 14 | non-small-molecule perturbagen (protein) |
| dtoxs_primary | Trastuzumab | antibody | 19 | non-small-molecule perturbagen (antibody) |
| temp_env | Blank | technical_blank | 32 | technical blank (concentration 0; neither control nor perturbation) |
| temp_env | ET-1 | peptide | 48 | non-small-molecule perturbagen (peptide) |

## View-level exclusions (in the corpus, filtered by the view tier)

<details><summary><code>dtoxs_cm.h5ad</code> (atc-moa): 1 compound entry (1 drug key) filtered by the tier</summary>

| compound | drug_key | reason |
|---|---|---|
| ZINC150345845 | `QMKSCTINHJSWPD-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |

</details>

- `screen_ctox.h5ad` (atc-moa): no view-level compound exclusions
<details><summary><code>temp_env_atc-moa.h5ad</code> (atc-moa): 348 compound entries (347 drug keys) filtered by the tier</summary>

| compound | drug_key | reason |
|---|---|---|
| 1,2,3-Benzotriazole | `QRUDEWIWKLJBPS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,2,3-Trichlorobenzene | `RELMFMZEBKVZJC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,2,4,5-Tetrachlorobenzene | `JHBKHLUZVFWLAG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,2-Benzisothiazolin-3-one | `DMSMPAJRVJJAGA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,2-Dinitrobenzene | `IZUKQUVSCNEFMJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,3-Diisopropylbenzene | `UNEATYXSUBPPKP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,3-Dinitrobenzene | `WDCYWAQPCXBPJA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,3-Diphenylguanidine | `OWRCNXZUPFZXOS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1,5,9-Cyclododecatriene | `ZOLLIQAKMYWTBR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 1-(Hydroxymethyl)-5,5-dimethylhydantoin | `SIQZJFKTROUNPI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,2',4,4'-Tetrahydroxybenzophenone | `WXNRYSGJLQFHBR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,2'-Bisphenol F | `MQCPOLNSJCWPGT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,2'-Methylenebis(4-methyl-6-tert-butylphenol) | `KGRVJHAUYBGFFP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4,4'-Trihydroxybenzophenone | `OKJFKPFBSPZTAH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4,5-Trichlorophenol | `LHJGJYXLEPZJPM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4,6-Trimethylphenol | `BPRYUXCVCCNUFE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4-Dichlorophenoxyacetic acid | `OVSKIKFHRZPJSS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4-Dihydroxybenzophenone | `ZXDDPOHVAMWLBH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4-Dinitrophenol | `UFBJCMHMOXMLKC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,4-Dinitrotoluene | `RMBFBMJGBANMMK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,5-Di-tert-butylbenzene-1,4-diol | `JZODKRWQWUWGCD-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,6,10-Trimethyl-2,6,10-triazaundecane | `SKCNNQDRNPQEFU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,6-Di-tert-butylphenol | `DKCPKDPYUFEZCP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,6-Diethylaniline | `FOYHNROGBXVLLX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,6-Dimethylaniline | `UFFBMTHBGFGIHF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,6-Dimethylphenol | `NXXYKOUNUYWIHA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2,6-Dinitrotoluene | `XTRDKALNCIHHNI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Anisidine | `VMPITZXILSNTON-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Methoxy-5-methylaniline | `WXWCDTXEKCVRRO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Methoxy-5-nitroaniline | `NIPDVSLAMPAWTP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Methyl-4,6-dinitrophenol | `ZXVONLUNISGICL-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Naphthylamine | `JBIJLHTVPXGSAM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Phenoxyethanol | `QCDWFXQBSFUVSP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Phenylphenol | `LLEMOWNGBBNAJR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-Tert-Butyl-5-methylphenol | `XOUQAVYLRNOXDO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-[(4-Hydroxyphenyl)methyl]phenol | `LVLNPXCISNPHLE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-[2-(2-Ethoxyethoxy)ethoxy]ethanol | `WFSMVVDJSNMRAR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 2-tert-Butyl-4-methoxyphenol | `MRBKEAMVRSLQPH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 3,3',5,5'-Tetrabromobisphenol A | `VEORPZCZECFIRK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 3,3'-Dimethylbenzidine | `NUIURNJTPRWVAP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 3,4-Dichloronitrobenzene | `NTBYINQTYWZXLH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 3-Hydroxyfluorene | `PVUBSZGNXLNTLX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 3-Nitrotoluene | `QZYHIOPPLUPUJF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4,4',4''-Ethane-1,1,1-triyltriphenol | `BRPSWMCDEYMRPE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4,4'-(9H-Fluorene-9,9-diyl)diphenol | `YWFPGFJLYRKYJZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4,4'-Diaminobiphenyl methane | `YBRVSVVVWCFQMG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4,4'-Methylenebis(N,N-dimethylaniline) | `JNRLEMMIVRBKJE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4,4'-Oxydianiline | `HLBLWEWZXPIGSM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4,4'-Sulfonyldiphenol | `VPWNQTHUCYMVMZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-(1,1,3,3-Tetramethylbutyl)phenol | `ISAVYTVYFVQUDY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-(2-Methylbutan-2-yl)phenol | `NRZWYNLTFLDQQX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-(Butan-2-yl)phenol | `ZUTYZAFDFLLILI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Androstene-3,17-dione | `AEMFNILZOJDQLW-QAGGRKNESA-N` | in corpus but outside the atc-moa tier |
| 4-Chloroaniline | `QSNSCYSYFYORTR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Cumylphenol | `QBDSZLJBMIMQRS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Hydroxy-3-methoxybenzaldehyde | `MWOOGOJBHIARFG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Nitrophenol | `BTJIUGUIPKRLHP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Nitrosodiphenylamine | `OIJHFHYPXWSVPF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Nitrotoluene | `ZPTVNYMJQHSSEA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-Octylphenol | `NTDQQZYCCIDJRK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 4-tert-Butylphenol | `QHPQWRBYOIRBIT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 5-Heptyldihydro-2(3H)-furanone | `PHXATPHONSXBIL-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 5HPP-33 | `LAKWINYVWJPHQW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 6:2 Fluorotelomer methacrylate | `CDXFIRXEAJABAZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| 9-Phenanthrol | `DZKIUEHLEXLYKM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acenaphthene | `CWRYPZZKDGJXCA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acenaphthylene | `HXGDTGSAIMULJN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acephate | `YASYVMFAVPKPKE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acetamiprid | `WCXDHFDTOYPNIE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acetochlor | `VTNQPKFIQCLBDU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acetophenone | `KWOLFJPFCHCOCG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acetyl tributyl citrate | `QZCLKYGREBVARF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Acibenzolar-S-methyl | `UELITFHSCLAHKR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Adipic acid | `WNLRTRBMVRJNCN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Aldrin | `QBYJBZPUGVGKQQ-SJJAEHHWSA-N` | in corpus but outside the atc-moa tier |
| Allethrin | `ZCVAOQKBXKSDMS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Allura Red C.I.16035 | `UQWIHFJXDRNUDP-UHFFFAOYSA-L` | in corpus but outside the atc-moa tier |
| Ametryn | `RQVYBGPQFYCBGX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Ammonium perfluorooctanoate | `SNGREZUHAYWORS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Anilazine | `IMHBYKMAHXWHRP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Anthracene | `MWPLVEDNUUSJAV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Asulam | `VGPYEHKOIGNJKV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Atrazine | `MXWJVTOOROXGIU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Azamethiphos | `VNKBTWQZTQIWDV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Azimilide | `MREBEPTUUMTTIA-PCLIKHOPSA-N` | in corpus but outside the atc-moa tier |
| Azinphos-methyl | `CJJOSEISRRTUQB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Azoxystrobin | `WFDXOXNFNRHQEC-GHRIWEEISA-N` | in corpus but outside the atc-moa tier |
| Benfluralin | `SMDHCQAYESWHAE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Benomyl | `RIOXQFHNBCKOKP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bentazone | `ZOMSMJKLGFBRBS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Benzidine | `HFACYLZERDEVSX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Benzophenone | `RWCCWEUUXYIKHB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Benzyl butyl phthalate | `IRIAEXORFWYRCZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bifenazate | `VHLKTXFWDRXILV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bifenthrin | `OMFRMAHOUUJSGP-IRHGGOMRSA-N` | in corpus but outside the atc-moa tier |
| Biphenyl | `ZUOUZKKEUPVFJK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bis(2-ethylhexyl) nonanedioate | `ZDWGXBPVPXVXMQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bis(2-ethylhexyl)hexanedioate | `SAOKZLXYCUGLFA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bisphenol A | `IISBACLAFKSPIT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bisphenol B | `HTVITOHKHWFJKO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Boscalid | `WYEMLYFITZORAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Bromacil | `CTSLUCNDVMMDHG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Buprofezin | `PRLVTUNWOQKEAI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Butylate | `BMTAFVWTTFSTOG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Butylparaben | `QFOHBWFCKVYLES-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| CJ-013610 | `VPTONMHDLLMOOV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| CP-122721 | `ZIWFCOIGUNPHPM-HKUYNNGSSA-N` | in corpus but outside the atc-moa tier |
| Candoxatril | `ZTWZVMIYIIVABD-OEMFJLHTSA-N` | in corpus but outside the atc-moa tier |
| Carbaryl | `CVXBEEMKQHEXEN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Carbendazim | `TWFZGCMQGLPBSX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Carbosulfan | `JLQUFIHWVLZVTJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Carboxin | `GYSSRZJIHXQEHQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Catechol | `YCIMNLLNPGFGHC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Chloridazon | `WYKYKTKDBLFHCY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Chloroneb | `PFIADAMVCJPXSF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Chlorpyrifos | `SBPBAQFWLVIOKP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Chlorpyrifos-methyl | `HRBKVYFZANMGRE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Chlorthal-dimethyl | `NPOJQCVWMSKXDN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Clodinafop-propargyl | `JBDHZKLJNAIJNC-LLVKDONJSA-N` | in corpus but outside the atc-moa tier |
| Clofentezine | `UXADOQPNKNTIHB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Clomazone | `KIEDNEWSYUYDSN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cloprop | `YNTJKQDWYXUTLZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Clorophene | `NCKMMSIFQUPKCK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Clothianidin | `PGOOBECODWQEAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Coumaphos | `BXNANOICGRISHX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cyanazine | `MZZBPDKVEFVLFF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cyazofamid | `YXKMMRDKEKCERS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cyclanilide | `GLWWLNJJJCTFMZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cycloate | `DFCAFRGABIXSDS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cyclopamine | `QASFUMOKHFSJGL-LAFRSMQTSA-N` | in corpus but outside the atc-moa tier |
| Cyproconazole | `UFNOUKDBUJZYDE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cyprodinil | `HAORKNGNJCEJBX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Cyromazine | `LVQDKIWDGQRHTE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Daminozide | `NOQGZXFMHARMLW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Darbufelone mesylate | `AKTXOQVMWSFEBQ-LCYFTJDESA-N` | in corpus but outside the atc-moa tier |
| Dazomet | `QAYICIQNSGETAS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Denatonium saccharide | `ZFQMTVNLDNXRNQ-UHFFFAOYSA-O` | in corpus but outside the atc-moa tier |
| Di(2-ethylhexyl) phthalate | `BJQHLKABXJIVAM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Di-n-octyl phthalate | `MQIUGAXCHLFZKX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diallyl phthalate | `QUDWYFHPNIMBFC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dibenzofuran | `TXCDCPKCNAJMEE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dibutyl adipate | `XTJFFFGAUHQWII-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dibutyl decanedioate | `PYGXAGIECVVIOZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dicamba | `IWEDIXLBFLAXBO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dichlobenil | `YOYAIZYFCNQIRF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dicloran | `BIXZHMJUSMUDOQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diclosulam | `QNXAVFXEJCPCJO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dicrotophos | `VEENJGZXVHKXNB-VOTSOKGWSA-N` | in corpus but outside the atc-moa tier |
| Dieldrin | `DFBKLUNHFCTMDC-PICURKEMSA-N` | in corpus but outside the atc-moa tier |
| Diethyl phthalate | `FLKPEMZONWLCSK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Difenoconazole | `BQYJATMQXGBDHF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Difenzoquat metilsulfate | `LBGPXIPGGRQBJW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diflubenzuron | `QQQYTWIFVNKMRW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diisobutyl phthalate | `MGWAVDBGNNKXQV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diisopropyl methylphosphonate | `WOAFDHWYKSOANX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dimethenamid | `JLYFCTQDENRSOL-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dimethoate | `MCWXGJITAZMZEV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dimethyl succinate | `MUXOBHXGJLMRAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diniconazole | `FBOUIAKEJMZPQG-AWNIVKPZSA-N` | in corpus but outside the atc-moa tier |
| Diphenolic acid | `VKOUCJUTMGHNOR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diphenylamine | `DMBHHRLKUKUOEG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dipropyl 2,5-pyridinedicarboxylate | `IITCWRFYJWUUPC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dithiopyr | `YUBJPYNSGLJZPQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Diuron | `XMTQQYYKAHVGBJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Dodecyltrimethylammonium chloride | `VICYBMUVWHJEFT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| EPN | `AIGRXSNSLVJMEA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Endothal | `GXEKYRXVRROBEV-FBXFSONDSA-N` | in corpus but outside the atc-moa tier |
| Endrin | `DFBKLUNHFCTMDC-GKRDHZSOSA-N` | in corpus but outside the atc-moa tier |
| Ethofumesate | `IRCMYGHHKLLGHV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Ethoprop | `VJYFKVYYMZPMAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Ethoxyquin | `DECIPOUIJURFOJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Etridiazole | `KQTVWCSONPJJPE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Eugenol | `RRAFCDWBNXTKKO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| FD&C Green No. 3 | `VVJKKKDJADMKNM-UHFFFAOYSA-L` | in corpus but outside the atc-moa tier |
| Famoxadone | `PCCSBWNGDMYFCW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenamiphos | `ZCJPOPBZHLUFHF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenarimol | `NHOWDZOIZKMVAI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenbuconazole | `RQDJADAKIFFEKQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenhexamid | `VDLGAVXLJYLFDH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenitrothion | `ZNOLGFHPUIJIMJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenoxaprop-ethyl | `PQKBPHSEKWERTG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenoxycarb | `HJUFTIJOISQSKQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenthion | `PNVJTZOFSHSLTO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fenvalerate | `NYPJDWWKZLNGGM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fipronil | `ZOCSXAVNDGMNBV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fluazifop-butyl | `VAIZTNZGPYBOGF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fluazinam | `UZCGKGPEKUCDTF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fludioxonil | `MUJOIMFVNIBMKC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Flufenacet | `IANUJLZYFUDJIH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Flumetralin | `PWNAWOCHVWERAR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Flumetsulam | `RXCPQSJAVKGONC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Flumioxazin | `FOUWCSDKDDHKQP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fluoranthene | `GVEPBJHOBDJJJI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fluorene | `NIHNNTQXNPWCJQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fluoxastrobin | `UFEODZBUAFNAEU-NLRVBDNBSA-N` | in corpus but outside the atc-moa tier |
| Flusilazole | `FQKUGOMFVDPBIZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Formetanate hydrochloride | `RMFNNCGOSPBBAD-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Fosthiazate | `DUFVKSUJRWYZQP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Genistein | `TZBJGXHYKVUXJN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Geranyl acetate | `HIGQPQRQIQDZMP-DHZHZOJOSA-N` | in corpus but outside the atc-moa tier |
| Glycidyl trimethylammonium chloride | `LRQFRDLKOBFOSG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Heptachlor | `FRCCEHPWNOQAEU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Hexachlorocyclopentadiene | `VUNCWTMEJYMOOR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Hexaconazole | `STMIIPIFODONDC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Hexythiazox | `XGWIJUOSCAQSSV-XHDPSFHLSA-N` | in corpus but outside the atc-moa tier |
| Hydramethylnon | `IQVNEKKDSLOHHK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Icaridin | `QLHULAHOXSSASE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Imazalil | `PZBPKYOVPCNPJY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Imazapic | `PVSGXWMWNRGTKE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Imazapyr | `CLQMBPJKHLGMQK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Imazethapyr | `XVOKUMIPKHGGTN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Imidacloprid | `YWTYJOPNNQFBPC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Iodosulfuron-methyl-sodium | `VWGAYSCWLXQJBQ-UHFFFAOYSA-M` | in corpus but outside the atc-moa tier |
| Iprodione | `ONUFESLQCSAYKA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Isoeugenol | `BJIOGJUNALELMI-ONEGZZNKSA-N` | in corpus but outside the atc-moa tier |
| Isopropalin | `NEKOXWSIMFDGMA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Isoxaben | `PMHURSZHKKJGBM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Isoxaflutole | `OYIKARCXOQLFHF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Lactofen | `CONWAEURSVPLRM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Lindane | `JLYXXMFPNIAWKQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Malaoxon | `WSORODGWGUUOBO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Maltol | `XPCTZQVDEJYUGT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Metalaxyl | `ZQEIXNIJLIKNTD-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methamidophos | `NNKVPIKMPCQWCG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methidathion | `MEBQXILRKZHVCX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methoxyfenozide | `QCAWEPFNJXQPAN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methyl (RS)-2-[4-(2,4-dichlorophenoxy)phenoxy]propionate | `BACHBFVBHLGWSL-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methyl 2-aminobenzoate | `VAMXMNNIEUEQDV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methyl dodecanoate | `UQDUPQYQJKYHQI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methyl parathion | `RLBIQVVOMOPOHC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methyleugenol | `ZYEMGPIYFIJGTP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Methylparaben | `LXCFILQKKLGQFO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Metolachlor | `WVQBLGZPHOPPFO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Mevinphos | `GEPDYQSQVLXLEU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Molinate | `DEDOPGXGGQYYMW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Monobenzyl phthalate | `XIKIUQUXDNHBFR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Monomethyl phthalate | `FNJSWIPFHMKRAT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Monuron | `BMLIZLVNXIYGCK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Myclobutanil | `HZJKXKUJVSEEFU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N,N'-Methylenebisacrylamide | `ZIUHHBKFKCYYJD-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N,N,4-Trimethylaniline | `GYVGXEWAOAAJEU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N,N-Diethylaniline | `GGSUCNLOZRCGPQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N-Ethyl-3-methylaniline | `GUYMMHOQXYZMJQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N-Ethylperfluorooctanesulfonamide | `CCEKAJIANROZEO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N-Nitrosodiphenylamine | `UBUCNCOMADRQHX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N-Phenyl-1,4-benzenediamine | `ATGUVEKSASEFFO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| N-acetyl-procainamide | `KEECCEWTUVWFCV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Naphthalene | `UFWIBTONFRDIAS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Napropamide | `WXZVAROIGSFCFJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Nitrapyrin | `DCUJJWWUNKIJPH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Nitrobenzene | `LQNUZADURLCDLV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Novaluron | `NJPPVKZQTLUDBO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Octabenzone | `QUAMTGJKVDWJEQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Octrizole | `IYAZLDLPUNDVAG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Octylbicycloheptenedicarboximide | `WLLGXSLBOPFWQV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Oryzalin | `UNAHYJYOSSSJHH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Oxadiazon | `CHNUNORXWHYHNE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Oxamyl | `KZAUOCCYDRDERY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Oxyfluorfen | `OQMBBFQZGJFLBU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Oxytetracycline hydrochloride | `OWFJMIVZYSDULZ-PXOLEDIWSA-N` | in corpus but outside the atc-moa tier |
| PHA-00568487 | `LUVXHMJTVXZFPD-ZDUSSCGKSA-N` | in corpus but outside the atc-moa tier |
| PK 11195 | `RAVIZVQZGXBOQO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Paclobutrazol | `RMOGWMIKYWRTKW-UONOGXRCSA-N` | in corpus but outside the atc-moa tier |
| Parathion | `LCCNCVORNKJIRZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pendimethalin | `CHIFOSRWCNZCFN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Penoxsulam | `SYJGKVOENHZYMQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pentachloronitrobenzene | `LKPLKUMXSAEKID-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pentachlorophenol | `IZUPBVBPLAPZRR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Perfluorodecanoic acid | `PCIUEQPBYFRTEM-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Perfluorohexanoic acid | `PXUULQAPEKKVAH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Perfluorononanoic acid | `UZUFPBIDKMEQEQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Perfluorooctanesulfonic acid | `YFSUTJLHUFNCNZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Perfluorooctanoic acid | `SNGREZUHAYWORS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Phenolphthalin | `FFFPYJTVNSSLBQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Phenothiazine | `WJFKNYWRSNBZNX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Phorate | `BULVZWIRKLYCBC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Phosalone | `IOUNQDKNJZEDEP-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Picloram | `NQQVFXUMIDALNH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pirimicarb | `YFGYUFNIOHWBOB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pirinixic acid | `SZRPDCCEHVWOJX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Potassium perfluorobutanesulfonate | `JGTNAGYHADQMCM-UHFFFAOYSA-M` | in corpus but outside the atc-moa tier |
| Potassium perfluorohexanesulfonate | `QZHDEAJFRJCDMF-UHFFFAOYSA-M` | in corpus but outside the atc-moa tier |
| Potassium perfluorooctanesulfonate | `YFSUTJLHUFNCNZ-UHFFFAOYSA-M` | in corpus but outside the atc-moa tier |
| Prochloraz | `TVLSRXXIMLFWEO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Prodiamine | `RSVPPPHXAASNOL-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Prometon | `ISEUFVQQFVOBCY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Propamocarb hydrochloride | `WZZLDXDUQPOXNW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Propanil | `LFULEKSKNZEWOE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Propetamphos | `BZNDWPRGXNILMS-VQHVLOKHSA-N` | in corpus but outside the atc-moa tier |
| Propiconazole | `STJLVHWMYQXCPB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Propoxycarbazone-sodium | `JTHMVYBOQLDDIY-UHFFFAOYSA-M` | in corpus but outside the atc-moa tier |
| Propyl gallate | `ZTHYODDOHIVTJV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Propylparaben | `QELSKZZBTMNZEB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Propyzamide | `PHNUZKMIPFFYSO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pymetrozine | `QHMTXANCGGJZRX-WUXMJOGZSA-N` | in corpus but outside the atc-moa tier |
| Pyraflufen-ethyl | `APTZNLHMIGJTEW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pyrene | `BBEAQIROQSPTKN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pyrimethanil | `ZLIBICFPKPWGIZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Pyriproxyfen | `NHDHVHZZCFYRSB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Quercetin | `REFJWTPEDVJJIY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Quinoline | `SMWDFEZZVXVKRB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Quinoxyfen | `WRPIRSINYZBGPK-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Resmethrin | `VEMKTZHHVJILDY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Rhodamine 6G | `IWWWBRIIGAXLCJ-BGABXYSRSA-N` | in corpus but outside the atc-moa tier |
| S-Bioallethrin | `ZCVAOQKBXKSDMS-PVAVHDDUSA-N` | in corpus but outside the atc-moa tier |
| Safrole | `ZMQAAUBTXCXRIC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Sematilide | `KHYPYQZQJSBPIX-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Sethoxydim | `CSPPKDPQLUUTND-NBVRZTHBSA-N` | in corpus but outside the atc-moa tier |
| Simazine | `ODCWYMIRDDJXKW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Spiroxamine | `PUYXTUJWRLOUCW-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Sucrose | `CZMRCDWAGMRECN-UGDNZRGBSA-N` | in corpus but outside the atc-moa tier |
| Sulisobenzone | `CXVGEDCSTKKODG-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tebuconazole | `PXMNMQRDXWABCY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tebuthiuron | `HBPDKDSFLXWOAE-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tefluthrin | `ZFHGXWPMULPQSE-SZGBIDFHSA-N` | in corpus but outside the atc-moa tier |
| Tepraloxydim | `IOYNQIMAUDJVEI-BMVIKAAMSA-N` | in corpus but outside the atc-moa tier |
| Terbacil | `NBQCNZYJJMBDKY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tetrabromophenolphthalein ethyl ester | `SQFXATUXPUCFFO-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Thiacloprid | `HOKKPVIRMVDYPB-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Thiamethoxam | `NWWZPOKUUAIXIW-DHZHZOJOSA-N` | in corpus but outside the atc-moa tier |
| Thiazopyr | `YIJZJEYQBAAWRJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Thidiazuron | `HFCYZXMHUIHAQI-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Thiobencarb | `QHTQREMOGMZHJV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Thiophanate-methyl | `QGHREAKMXXNCOA-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tralkoxydim | `DQFPEYARZIQXRM-LTGZKZEYSA-N` | in corpus but outside the atc-moa tier |
| Tri-allate | `MWBPRDONLNQCFV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triadimefon | `WURBVZBTWMNKQT-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triadimenol | `BAZVSMNPJJMILC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tribufos | `ZOKXUAHZSKEQSS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tributyl phosphate | `STCOOQWBFONSKY-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triclopyr | `REEQLXCGVXDJSQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triethylene glycol bis(2-ethylhexanoate) | `FRQDZJMEHSJOPU-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Trifloxystrobin | `ONCZDRURRATYFI-TVJDWZFNSA-N` | in corpus but outside the atc-moa tier |
| Trifloxysulfuron-sodium | `AIMMSOZBPYFASU-UHFFFAOYSA-M` | in corpus but outside the atc-moa tier |
| Triflumizole | `HSMVPDGQOIQYSR-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triglycidyl isocyanurate | `OUPZKGBUJRBPGC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Trinexapac-ethyl | `RVKCCVTVZORVGD-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triphenyl phosphate | `XZZNDPSIHUTMOC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tris(1,3-dichloro-2-propyl) phosphate | `ASLWPAWFJZFCKF-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Tris(2-chloroethyl) phosphate | `HQUQLFOMPYWACS-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| Triticonazole | `PPDBOQMNKNNODG-NTEUORMPSA-N` | in corpus but outside the atc-moa tier |
| Vernolate | `CCPPLLJZDQAOHD-FLIBITNWSA-M` | in corpus but outside the atc-moa tier |
| alpha-Isomethylionone | `JRJBVWJSTHECJK-PKNBQFBNSA-N` | in corpus but outside the atc-moa tier |
| o-Cresol | `QWVGKYWNOKOFNN-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| p,p'-DDD | `AHJKRLASYNVKDZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| p-Cresol | `IWDCLRJOBJJRNH-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| p-Xylenol blue | `MGUKYHHAGPFJMC-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |
| tert-Butylhydro-quinone (TBHQ) | `BGNXCDMCOKJUMV-UHFFFAOYSA-N` | in corpus but outside the atc-moa tier |

</details>

<details><summary><code>temp_env_atc-moa-mesh.h5ad</code> (atc-moa-mesh): 209 compound entries (208 drug keys) filtered by the tier</summary>

| compound | drug_key | reason |
|---|---|---|
| 1,2,3-Benzotriazole | `QRUDEWIWKLJBPS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,2,3-Trichlorobenzene | `RELMFMZEBKVZJC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,2,4,5-Tetrachlorobenzene | `JHBKHLUZVFWLAG-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,2-Dinitrobenzene | `IZUKQUVSCNEFMJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,3-Diisopropylbenzene | `UNEATYXSUBPPKP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,3-Dinitrobenzene | `WDCYWAQPCXBPJA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,3-Diphenylguanidine | `OWRCNXZUPFZXOS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1,5,9-Cyclododecatriene | `ZOLLIQAKMYWTBR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 1-(Hydroxymethyl)-5,5-dimethylhydantoin | `SIQZJFKTROUNPI-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,2',4,4'-Tetrahydroxybenzophenone | `WXNRYSGJLQFHBR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,2'-Bisphenol F | `MQCPOLNSJCWPGT-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,2'-Methylenebis(4-methyl-6-tert-butylphenol) | `KGRVJHAUYBGFFP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,4,4'-Trihydroxybenzophenone | `OKJFKPFBSPZTAH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,4,5-Trichlorophenol | `LHJGJYXLEPZJPM-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,4,6-Trimethylphenol | `BPRYUXCVCCNUFE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,4-Dihydroxybenzophenone | `ZXDDPOHVAMWLBH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,6,10-Trimethyl-2,6,10-triazaundecane | `SKCNNQDRNPQEFU-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,6-Di-tert-butylphenol | `DKCPKDPYUFEZCP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,6-Diethylaniline | `FOYHNROGBXVLLX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,6-Dimethylaniline | `UFFBMTHBGFGIHF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,6-Dimethylphenol | `NXXYKOUNUYWIHA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2,6-Dinitrotoluene | `XTRDKALNCIHHNI-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2-Anisidine | `VMPITZXILSNTON-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2-Methoxy-5-methylaniline | `WXWCDTXEKCVRRO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2-Methoxy-5-nitroaniline | `NIPDVSLAMPAWTP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2-Tert-Butyl-5-methylphenol | `XOUQAVYLRNOXDO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2-[(4-Hydroxyphenyl)methyl]phenol | `LVLNPXCISNPHLE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 2-[2-(2-Ethoxyethoxy)ethoxy]ethanol | `WFSMVVDJSNMRAR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 3,3',5,5'-Tetrabromobisphenol A | `VEORPZCZECFIRK-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 3,3'-Dimethylbenzidine | `NUIURNJTPRWVAP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 3,4-Dichloronitrobenzene | `NTBYINQTYWZXLH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 3-Hydroxyfluorene | `PVUBSZGNXLNTLX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 3-Nitrotoluene | `QZYHIOPPLUPUJF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4,4',4''-Ethane-1,1,1-triyltriphenol | `BRPSWMCDEYMRPE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4,4'-(9H-Fluorene-9,9-diyl)diphenol | `YWFPGFJLYRKYJZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4,4'-Methylenebis(N,N-dimethylaniline) | `JNRLEMMIVRBKJE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4,4'-Oxydianiline | `HLBLWEWZXPIGSM-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4,4'-Sulfonyldiphenol | `VPWNQTHUCYMVMZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-(2-Methylbutan-2-yl)phenol | `NRZWYNLTFLDQQX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-(Butan-2-yl)phenol | `ZUTYZAFDFLLILI-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-Androstene-3,17-dione | `AEMFNILZOJDQLW-QAGGRKNESA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-Chloroaniline | `QSNSCYSYFYORTR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-Cumylphenol | `QBDSZLJBMIMQRS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-Nitrophenol | `BTJIUGUIPKRLHP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-Nitrosodiphenylamine | `OIJHFHYPXWSVPF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-Nitrotoluene | `ZPTVNYMJQHSSEA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 4-tert-Butylphenol | `QHPQWRBYOIRBIT-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 5-Heptyldihydro-2(3H)-furanone | `PHXATPHONSXBIL-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 5HPP-33 | `LAKWINYVWJPHQW-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| 6:2 Fluorotelomer methacrylate | `CDXFIRXEAJABAZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Acenaphthene | `CWRYPZZKDGJXCA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Acenaphthylene | `HXGDTGSAIMULJN-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Acetamiprid | `WCXDHFDTOYPNIE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Acetyl tributyl citrate | `QZCLKYGREBVARF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Acibenzolar-S-methyl | `UELITFHSCLAHKR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Adipic acid | `WNLRTRBMVRJNCN-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Allura Red C.I.16035 | `UQWIHFJXDRNUDP-UHFFFAOYSA-L` | in corpus but outside the atc-moa-mesh tier |
| Ammonium perfluorooctanoate | `SNGREZUHAYWORS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Anilazine | `IMHBYKMAHXWHRP-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Anthracene | `MWPLVEDNUUSJAV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Asulam | `VGPYEHKOIGNJKV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Azoxystrobin | `WFDXOXNFNRHQEC-GHRIWEEISA-N` | in corpus but outside the atc-moa-mesh tier |
| Benfluralin | `SMDHCQAYESWHAE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Benzidine | `HFACYLZERDEVSX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Bifenazate | `VHLKTXFWDRXILV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Bifenthrin | `OMFRMAHOUUJSGP-IRHGGOMRSA-N` | in corpus but outside the atc-moa-mesh tier |
| Bis(2-ethylhexyl) nonanedioate | `ZDWGXBPVPXVXMQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Bis(2-ethylhexyl)hexanedioate | `SAOKZLXYCUGLFA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Bisphenol B | `HTVITOHKHWFJKO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Boscalid | `WYEMLYFITZORAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Buprofezin | `PRLVTUNWOQKEAI-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Butylate | `BMTAFVWTTFSTOG-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Butylparaben | `QFOHBWFCKVYLES-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| CJ-013610 | `VPTONMHDLLMOOV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| CP-122721 | `ZIWFCOIGUNPHPM-HKUYNNGSSA-N` | in corpus but outside the atc-moa-mesh tier |
| Candoxatril | `ZTWZVMIYIIVABD-OEMFJLHTSA-N` | in corpus but outside the atc-moa-mesh tier |
| Carbosulfan | `JLQUFIHWVLZVTJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Catechol | `YCIMNLLNPGFGHC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Chloridazon | `WYKYKTKDBLFHCY-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Chloroneb | `PFIADAMVCJPXSF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Chlorpyrifos-methyl | `HRBKVYFZANMGRE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Chlorthal-dimethyl | `NPOJQCVWMSKXDN-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Clodinafop-propargyl | `JBDHZKLJNAIJNC-LLVKDONJSA-N` | in corpus but outside the atc-moa-mesh tier |
| Clofentezine | `UXADOQPNKNTIHB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Clomazone | `KIEDNEWSYUYDSN-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Cloprop | `YNTJKQDWYXUTLZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Cyazofamid | `YXKMMRDKEKCERS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Cyclanilide | `GLWWLNJJJCTFMZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Cycloate | `DFCAFRGABIXSDS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Cyclopamine | `QASFUMOKHFSJGL-LAFRSMQTSA-N` | in corpus but outside the atc-moa-mesh tier |
| Cyprodinil | `HAORKNGNJCEJBX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Darbufelone mesylate | `AKTXOQVMWSFEBQ-LCYFTJDESA-N` | in corpus but outside the atc-moa-mesh tier |
| Dazomet | `QAYICIQNSGETAS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Denatonium saccharide | `ZFQMTVNLDNXRNQ-UHFFFAOYSA-O` | in corpus but outside the atc-moa-mesh tier |
| Di(2-ethylhexyl) phthalate | `BJQHLKABXJIVAM-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Di-n-octyl phthalate | `MQIUGAXCHLFZKX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diallyl phthalate | `QUDWYFHPNIMBFC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dibenzofuran | `TXCDCPKCNAJMEE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dibutyl adipate | `XTJFFFGAUHQWII-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dibutyl decanedioate | `PYGXAGIECVVIOZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dicloran | `BIXZHMJUSMUDOQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diclosulam | `QNXAVFXEJCPCJO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dicrotophos | `VEENJGZXVHKXNB-VOTSOKGWSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dieldrin | `DFBKLUNHFCTMDC-PICURKEMSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diethyl phthalate | `FLKPEMZONWLCSK-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Difenzoquat metilsulfate | `LBGPXIPGGRQBJW-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diflubenzuron | `QQQYTWIFVNKMRW-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diisobutyl phthalate | `MGWAVDBGNNKXQV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diisopropyl methylphosphonate | `WOAFDHWYKSOANX-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dimethyl succinate | `MUXOBHXGJLMRAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diniconazole | `FBOUIAKEJMZPQG-AWNIVKPZSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diphenolic acid | `VKOUCJUTMGHNOR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Diphenylamine | `DMBHHRLKUKUOEG-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dipropyl 2,5-pyridinedicarboxylate | `IITCWRFYJWUUPC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Dithiopyr | `YUBJPYNSGLJZPQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| EPN | `AIGRXSNSLVJMEA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Endothal | `GXEKYRXVRROBEV-FBXFSONDSA-N` | in corpus but outside the atc-moa-mesh tier |
| Ethoprop | `VJYFKVYYMZPMAB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Ethoxyquin | `DECIPOUIJURFOJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Etridiazole | `KQTVWCSONPJJPE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| FD&C Green No. 3 | `VVJKKKDJADMKNM-UHFFFAOYSA-L` | in corpus but outside the atc-moa-mesh tier |
| Famoxadone | `PCCSBWNGDMYFCW-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Fenhexamid | `VDLGAVXLJYLFDH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Fluazinam | `UZCGKGPEKUCDTF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Fludioxonil | `MUJOIMFVNIBMKC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Flufenacet | `IANUJLZYFUDJIH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Flumetralin | `PWNAWOCHVWERAR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Flumetsulam | `RXCPQSJAVKGONC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Formetanate hydrochloride | `RMFNNCGOSPBBAD-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Geranyl acetate | `HIGQPQRQIQDZMP-DHZHZOJOSA-N` | in corpus but outside the atc-moa-mesh tier |
| Glycidyl trimethylammonium chloride | `LRQFRDLKOBFOSG-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Hexaconazole | `STMIIPIFODONDC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Hexythiazox | `XGWIJUOSCAQSSV-XHDPSFHLSA-N` | in corpus but outside the atc-moa-mesh tier |
| Iodosulfuron-methyl-sodium | `VWGAYSCWLXQJBQ-UHFFFAOYSA-M` | in corpus but outside the atc-moa-mesh tier |
| Isoeugenol | `BJIOGJUNALELMI-ONEGZZNKSA-N` | in corpus but outside the atc-moa-mesh tier |
| Isopropalin | `NEKOXWSIMFDGMA-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Isoxaben | `PMHURSZHKKJGBM-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Isoxaflutole | `OYIKARCXOQLFHF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Maltol | `XPCTZQVDEJYUGT-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Methyl (RS)-2-[4-(2,4-dichlorophenoxy)phenoxy]propionate | `BACHBFVBHLGWSL-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Methyl 2-aminobenzoate | `VAMXMNNIEUEQDV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Methyl dodecanoate | `UQDUPQYQJKYHQI-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Monobenzyl phthalate | `XIKIUQUXDNHBFR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Monomethyl phthalate | `FNJSWIPFHMKRAT-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Monuron | `BMLIZLVNXIYGCK-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| N,N'-Methylenebisacrylamide | `ZIUHHBKFKCYYJD-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| N,N,4-Trimethylaniline | `GYVGXEWAOAAJEU-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| N,N-Diethylaniline | `GGSUCNLOZRCGPQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| N-Ethyl-3-methylaniline | `GUYMMHOQXYZMJQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| N-Ethylperfluorooctanesulfonamide | `CCEKAJIANROZEO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| N-Phenyl-1,4-benzenediamine | `ATGUVEKSASEFFO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Naphthalene | `UFWIBTONFRDIAS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Napropamide | `WXZVAROIGSFCFJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Nitrapyrin | `DCUJJWWUNKIJPH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Nitrobenzene | `LQNUZADURLCDLV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Octabenzone | `QUAMTGJKVDWJEQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Octrizole | `IYAZLDLPUNDVAG-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Octylbicycloheptenedicarboximide | `WLLGXSLBOPFWQV-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Oxadiazon | `CHNUNORXWHYHNE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Oxamyl | `KZAUOCCYDRDERY-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Oxyfluorfen | `OQMBBFQZGJFLBU-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Oxytetracycline hydrochloride | `OWFJMIVZYSDULZ-PXOLEDIWSA-N` | in corpus but outside the atc-moa-mesh tier |
| PHA-00568487 | `LUVXHMJTVXZFPD-ZDUSSCGKSA-N` | in corpus but outside the atc-moa-mesh tier |
| Paclobutrazol | `RMOGWMIKYWRTKW-UONOGXRCSA-N` | in corpus but outside the atc-moa-mesh tier |
| Penoxsulam | `SYJGKVOENHZYMQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Perfluorodecanoic acid | `PCIUEQPBYFRTEM-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Perfluorohexanoic acid | `PXUULQAPEKKVAH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Perfluorooctanesulfonic acid | `YFSUTJLHUFNCNZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Perfluorooctanoic acid | `SNGREZUHAYWORS-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Phenolphthalin | `FFFPYJTVNSSLBQ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Potassium perfluorobutanesulfonate | `JGTNAGYHADQMCM-UHFFFAOYSA-M` | in corpus but outside the atc-moa-mesh tier |
| Potassium perfluorohexanesulfonate | `QZHDEAJFRJCDMF-UHFFFAOYSA-M` | in corpus but outside the atc-moa-mesh tier |
| Potassium perfluorooctanesulfonate | `YFSUTJLHUFNCNZ-UHFFFAOYSA-M` | in corpus but outside the atc-moa-mesh tier |
| Prodiamine | `RSVPPPHXAASNOL-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Prometon | `ISEUFVQQFVOBCY-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Propamocarb hydrochloride | `WZZLDXDUQPOXNW-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Propetamphos | `BZNDWPRGXNILMS-VQHVLOKHSA-N` | in corpus but outside the atc-moa-mesh tier |
| Propiconazole | `STJLVHWMYQXCPB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Propyzamide | `PHNUZKMIPFFYSO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Pymetrozine | `QHMTXANCGGJZRX-WUXMJOGZSA-N` | in corpus but outside the atc-moa-mesh tier |
| Pyraflufen-ethyl | `APTZNLHMIGJTEW-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Pyrimethanil | `ZLIBICFPKPWGIZ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Quinoline | `SMWDFEZZVXVKRB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Quinoxyfen | `WRPIRSINYZBGPK-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Rhodamine 6G | `IWWWBRIIGAXLCJ-BGABXYSRSA-N` | in corpus but outside the atc-moa-mesh tier |
| Safrole | `ZMQAAUBTXCXRIC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Sethoxydim | `CSPPKDPQLUUTND-NBVRZTHBSA-N` | in corpus but outside the atc-moa-mesh tier |
| Sulisobenzone | `CXVGEDCSTKKODG-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Tebuthiuron | `HBPDKDSFLXWOAE-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Tefluthrin | `ZFHGXWPMULPQSE-SZGBIDFHSA-N` | in corpus but outside the atc-moa-mesh tier |
| Tepraloxydim | `IOYNQIMAUDJVEI-BMVIKAAMSA-N` | in corpus but outside the atc-moa-mesh tier |
| Terbacil | `NBQCNZYJJMBDKY-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Tetrabromophenolphthalein ethyl ester | `SQFXATUXPUCFFO-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Thiacloprid | `HOKKPVIRMVDYPB-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Thiazopyr | `YIJZJEYQBAAWRJ-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Thidiazuron | `HFCYZXMHUIHAQI-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Tralkoxydim | `DQFPEYARZIQXRM-LTGZKZEYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Triethylene glycol bis(2-ethylhexanoate) | `FRQDZJMEHSJOPU-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Trifloxysulfuron-sodium | `AIMMSOZBPYFASU-UHFFFAOYSA-M` | in corpus but outside the atc-moa-mesh tier |
| Triflumizole | `HSMVPDGQOIQYSR-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Trinexapac-ethyl | `RVKCCVTVZORVGD-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Triphenyl phosphate | `XZZNDPSIHUTMOC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Tris(1,3-dichloro-2-propyl) phosphate | `ASLWPAWFJZFCKF-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| Triticonazole | `PPDBOQMNKNNODG-NTEUORMPSA-N` | in corpus but outside the atc-moa-mesh tier |
| Vernolate | `CCPPLLJZDQAOHD-FLIBITNWSA-M` | in corpus but outside the atc-moa-mesh tier |
| alpha-Isomethylionone | `JRJBVWJSTHECJK-PKNBQFBNSA-N` | in corpus but outside the atc-moa-mesh tier |
| o-Cresol | `QWVGKYWNOKOFNN-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| p-Cresol | `IWDCLRJOBJJRNH-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |
| p-Xylenol blue | `MGUKYHHAGPFJMC-UHFFFAOYSA-N` | in corpus but outside the atc-moa-mesh tier |

</details>

- `temp_env_all-smiles.h5ad` (all-smiles): no view-level compound exclusions

## Salt / free-base collapses

These names share one `drug_key` (parent structure) but remain separate
rows, so a cold-drug split cannot leak across them:

- Tamoxifen / Tamoxifen citrate: `NKANXQFJJICGDU-QPLCGJKRSA-N`
- Diltiazem / Diltiazem hydrochloride: `HSUGRBWQSSZJOP-RTWAWAEBSA-N`
- Estradiol / 17beta-Estradiol: `VOXZDWNPVJITMN-ZBRFXRBCSA-N`
- Perfluorooctanoic acid / Ammonium perfluorooctanoate: `SNGREZUHAYWORS-UHFFFAOYSA-N`

Doxorubicin and epirubicin are 4'-epimers and carry **different** parent
keys (connectivity-only keys would collapse them; the build rejects
14-character keys).

## Tier counts

| tier | drug keys (corpus) |
|---|---|
| atc-moa | 177 |
| atc-moa-mesh | 316 |
| all-smiles | 525 |

## Drift

Rerun matched `corpus/v1/` exactly: no gene or drug drift in any source.

## Manifest

| view | obs | perturbed | controls | drug keys | genes | content sha256 |
|---|---|---|---|---|---|---|
| `dtoxs_cm.h5ad` | 1,075 | 957 | 118 | 46 | 18,164 | `17198963e925c742` |
| `screen_ctox.h5ad` | 3,072 | 2,520 | 552 | 42 | 18,164 | `c0fffab49d0a2a91` |
| `temp_env_atc-moa.h5ad` | 1,319 | 1,047 | 272 | 113 | 18,164 | `91591d58dd63aed0` |
| `temp_env_atc-moa-mesh.h5ad` | 2,582 | 2,310 | 272 | 252 | 18,164 | `4633eaef47ea3860` |
| `temp_env_all-smiles.h5ad` | 4,478 | 4,206 | 272 | 460 | 18,164 | `f84066c4dfd2c017` |

## Caveats

- `dtoxs_primary` doses are literature/clinical Cmax values; `DIC`
  (diclofenac, 3 wells) uses the 5 uM recorded in the DToxS DCIC design
  workbook and is marked `dose_source=measured`.
- `screen_ctox` doses are derived from the curated paper Table-1 top
  concentration (`dose_source=derived_top_dose`), not measured.
- The DToxS `AMI` arm is a pentyl-benzofuran homolog (CID 16218975), not
  amiodarone; it is in `all-smiles` only.
- `dtoxs_cm.h5ad` intentionally mixes `primary_cm` and `ipsc_cm`; cell
  type and donor are covariates, not collapsed constants.
- The 32 `temp_env` `Blank` wells are technical blanks, neither
  controls nor perturbations.
