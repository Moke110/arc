# LINCS2020 `processed_data.h5` 审计报告

## 文件与结构

- **路径**: `training/data/LINCS2020_set/processed_data.h5`
- **规模**: 90184 条样本 × 978 基因（`x1`, `x2` 为 `float32`）
- **其它字段**: `canonical_smiles`, `cid`, `sig`（object）

## 关于「NaN」的结论

对 **全量** `x1`、`x2` 扫描结果：

| 检查项 | x1 | x2 |
|--------|----|----|
| `NaN` 元素个数 | **0** | **0** |
| `Inf` 元素个数 | **0** | **0** |
| 含任意 `NaN` 的行数 | **0** | **0** |

**结论：HDF5 中不存在缺失值（NaN / Inf）。** 训练日志里 `x1_rec_pearson` 等为 `nan` 的成因**不是**原始数据里的 NaN。

## 真正的问题：常数谱（基因维方差为 0）

`scipy.stats.pearsonr` 在任一向量沿基因维标准差为 **0** 时，Pearson **未定义**，结果为 `nan`。  
若按样本累加后求平均，**任意一个样本** 出现 `nan` 会把整体均值污染为 `nan`。

本数据集中存在：

- **`x1` 方差为 0 的行数**: **6**（行索引 `27505, 27515, 27689, 27703, 28823, 28828`）
- **`x2` 方差为 0 的行数**: **3**（与 `x1` 中 3 行重叠：JURKAT 上全 0 等）

典型情况：JURKAT 上 `x1`/`x2` 为常数 **0**；THP1 上部分为常数 **≈4.7667**（整条 profile 相同）。

## 与 split 的交叉

| 划分 | 命中「x1 或 x2 为零方差」的 LINCS 行 |
|------|--------------------------------------|
| train | 5 |
| valid | 1 |
| **test** | **0** |

**说明**：测试集 **不** 引用上述零方差行；若仅出现 **测试集** 上全谱 Pearson 为 `nan`，则更可能来自 **模型输出**（`x1_rec` 等）在部分样本上为常数或数值退化，与 `pearson` 聚合方式有关。  
对 **训练/验证** 中仍引用这些行的样本，零方差谱会直接导致相关 Pearson 未定义。

## 已做修复（本仓库）

1. **`src/utils.py` 中 `pearson()`**  
   在任一侧标准差为 0 时给出**可聚合**的约定：两向量近似相等则 `1.0`，否则 `0.0`；避免整条测试指标被 `nan` 污染。

2. **本目录脚本（可选数据层修复）**  
   - `audit_lincs2020_h5.py`：重复上述审计并写出 `LINCS2020_bad_row_indices.json`  
   - `fix_lincs2020_constant_profiles.py`：对零方差行加极小高斯抖动（默认 `eps=1e-7`），可生成 `processed_data_fixed.h5`

### 生成修复版 HDF5（约 700MB）

```bash
python3 .cursor/fix_lincs2020_constant_profiles.py -o .cursor/processed_data_fixed.h5
```

原地替换（先备份 `.bak`）：

```bash
python3 .cursor/fix_lincs2020_constant_profiles.py \
  --input training/data/LINCS2020_set/processed_data.h5 \
  --output training/data/LINCS2020_set/processed_data.h5 \
  --in-place
```

## 产出文件清单（`.cursor/`）

| 文件 | 说明 |
|------|------|
| `audit_lincs2020_h5.py` | 审计脚本 |
| `fix_lincs2020_constant_profiles.py` | 零方差行抖动修复脚本 |
| `LINCS2020_bad_row_indices.json` | 问题行索引（由审计脚本生成） |
| `LINCS2020_processed_data_audit_report.md` | 本报告 |

