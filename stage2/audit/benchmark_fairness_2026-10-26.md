# Stage 2 benchmark fairness audit — 2026-10-26

## 结论与审计边界

现有结果支持“PRnet 是当前 5k 配方下值得优先推进的候选”，不支持“已公平证明 PRnet 在各模型充分发挥后最好”。统一打分基本可靠；模型输入/输出语义、正常收敛、原始配方保真度、调优充分性仍存在实质缺口。最优先的问题是 **5k PerturbNet 的归一化分母不一致**。

本次审计读取了本地部署记录、ADR 0027/0028、数据构建代码、七个适配器相关代码、保真度记录、训练记录及保存的预测。未重训模型，未修改历史榜单或冻结协议；也未独立复现每篇论文的最佳结果。下面的作者配方比较依据本地 vendored code 和已有 fidelity ledger，不代表重新核验了所有论文实验。

## 已独立验证的部分

- 18k 的现有测试 46/46 通过；5k 为 53/53（2026-10-26）。包含数据哈希、分组、划分、Δ、测试 bundle 无表达字段、指标 sanity 等；其中部分适配器测试属于配置/静态契约检查，不能替代真实模型语义验收。
- 对两套实验的 PRnet、PerturbNet、delta_mean、Squidiff 共 40 个 fold，从 `pred_test.npy` 和组 ID 重新聚合，独立构造 top/bottom-200 标签并调用 sklearn AP；与 leaderboard 主指标最大差异小于 1e-15。
- 两套实验共享 well、drug split 与组定义，预测先按 well 汇总至 drug-dose 组，再作方向宏平均、组宏平均、fold 平均。各药均有三个剂量，因此这里组宏平均与先药内平均再药间平均相容。
- test bundle 不携带 treated expression；CPA 的 test 表达占位为零，注册测试药物结构不等于泄露其扰动表达。未在检查到的适配路径发现直接读取测试 treated truth 来拟合的行为。
- 实际 basal bundle 是 16 × 12 = **192 个 DMSO wells**。源数据的 272 controls 还包含其他 control 类型，不应把“272 controls”写成“272 DMSO”。同板对照共享是协议允许的信息，不等于测试药物表达泄漏。

| 实验 | PRnet AP@200 | PerturbNet AP@200 | delta_mean AP@200 | 含义 |
|---|---:|---:|---:|---|
| 18,164 genes | 0.1643265 | 0.1860556 | 0.1988713 | 没有模型超过均值基线；PerturbNet 为最佳模型行 |
| 5,000 genes | 0.2734561 | 0.1928778 | 0.2491443 | PRnet 在现有配方下领先；PerturbNet 行有尺度缺陷 |

证据：`../03_models_cmbench/03.{1_cold_drug_18k,2_hvg_sensitivity_5k}/execution/script/{metrics,evaluate,cases,build_data}.py`、各自 `README.md`、`output/results/leaderboard.csv`、原始 `output/units/`。

## 1. 已确认：5k PerturbNet 预测与评测真值不在同一归一化空间

5k `scripts/build_data.py:274-284` 先验证全 18,164 基因的 cp10k，再截取所选列。因此 5k 真值保留的是 **全基因库大小为分母** 的表达。

PerturbNet `scripts/adapters/perturbnet/run.py:339` 用 5k counts 训练 scVI；`:440` 的库大小也仅为这 5k counts 之和；`:707-711` 将 decoder rate 除以该库大小再乘 10,000，然后取 log1p、减去全基因分母下的 control mean。产生的是 **5k 内重新闭合至 10,000** 的表达，无法与冻结真值直接相减。

实测保存的五折 test 预测，加回同板 control mean、expm1 后：

| 对象 | 5k cp10k 总量范围 | 中位数 |
|---|---:|---:|
| 所有 4,478 个 source wells 的冻结真值 | 2,211.45–5,281.74 | 2,918.42 |
| PerturbNet 五折 4,206 个 test wells 的预测 | 9,999.995–10,000.006 | 9,999.999 |

所选基因占全基因 counts 的比例为 22.11%–52.82%，中位数 29.18%。该问题不是“少量浮点误差”，也不是换一个 median library 即可解决：现有代码先乘后除同一 library，尺度会抵消。

影响：RMSE/R² 明显受影响；AUPRC 也**不保证不变**，因为在绝对表达 rescale 后再 log1p 并减去逐基因不同的 control，不是对最终 Δ 做共同的单调变换。现有结果不能区分模型本身变化与这个错误的贡献，也不能预先断言修复后会超过 PRnet。

建议：新版本必须明确保留全基因分母；可以评估全轴 decoder 后切片的参考路径，或设计能够建模所选基因相对全轴质量的输出桥接。不要用 test treated library/质量比例进行事后校准。重算/重训受影响单元并保留旧结果作为 provenance。

另一个应验收的估计量问题：代码先对 8 个 flow draws 的线性表达取平均再 log1p，而 benchmark target 是各 well 的 log1p 后再平均，两者一般不相等；ledger 将 ZINB 的 `px_rate` 直接称为 predictive mean 也需要检查 zero-inflation 语义。此处不据此量化额外损失。

## 2. 已声明但不能视为“实力已发挥”的配方限制

| 模型 | 当前关键限制 | 可以下的结论 |
|---|---|---|
| PRnet | 5k 宽度贴近已有配方；从零训练、小样本；以 DMSO 均值为输入、零噪声点预测 | 当前 5k 排名可靠，但非全基因或疗效优胜证明 |
| CPA | 使用结构 embedding 冷药路径；未见药物 dose scaler 无该药监督；原生 `cpa_metric` 选 checkpoint | 是该扩展配置的结果，仍需验证 cold counterfactual 校准与选择目标 |
| chemCPA | 省略 L1000 pretraining，5k from scratch；五折最佳 epoch 为 201/201/199/200/201 | 不能据此否定完整迁移学习流程；端点最优提示预算充分性未证实 |
| biolord | 从已注册 DMSO latent 作 counterfactual；每 50 epoch 一次选点；四折最佳点在 500 上限 | 可以评估当前共享对照情境，不能宣称新 basal context 的 inductive 泛化；收敛/预算需逐折看 |
| PerturbNet | ZINC 大规模 ChemicalVAE 训练被约 220 个 train SMILES 替代；部分 SMILES 不可编码；不训练 untreated counts | 从零小样本配方受到强约束，另叠加上述 5k 输出错误 |
| Squidiff | 作者 Sci-Plex 药物结构示例使用 200 genes，此处为该示例的 25 倍；100k steps 后预测 RMSE 约 74 | 未通过正常工作验收，不能当成完整能力的低排名；200 不是架构上限 |
| CMonge | 1k→5k，50-d latent；all-drug 500k iteration 配方改为 10k；每条件约 3 wells | 是受资源/数据限制的适配结果，不能称充分训练的分布建模比较 |

原始来源：`03.2_hvg_sensitivity_5k/execution/models/<model>/FIDELITY.md` 与 `output/results/run_cards.csv`。

Squidiff 补充核查：作者 [README](https://github.com/siyuh/Squidiff#training-squidiff) 同时提供 500-gene 通用训练示例与 200-gene 药物结构示例。本地 `pert_model/Squidiff/source/Squidiff/MLPModel.py:134-135` 按 `gene_size` / `output_dim` 动态构造输入输出层，没有 200-gene 硬限制。因此只能说本次 5k 偏离了该药物示例，不能说该模型设计上只能处理 200 genes；已有 ledger 中“never validated ... above 200”的绝对表述也未获本次核查支持。当前不收敛的根因尚未确定，不能归因于越过所谓 200-gene 上限。

### Squidiff 全转录组资源复核

用户明确最终输出要求为 15,000+ genes；不能把仅在小 HVG 面板成功作为最终交付。2026-10-26 直接在 torch meta device 上实例化 vendored `MLPModel`，保持当前三层、2048 hidden、drug structure encoder 配置，仅改变输入输出基因数，计数如下（无需分配实际权重或重训）：

| genes | 参数量 | 相对 200 genes | FP32 参数本体 MB |
|---:|---:|---:|---:|
| 200 | 54,221,060 | 1.00× | 216.88 |
| 5,000 | 83,717,060 | 1.54× | 334.87 |
| 15,000 | 145,167,060 | 2.68× | 580.67 |
| 18,164 | 164,609,840 | 3.04× | 658.44 |

固定 hidden 时每增加一个基因增加 6,145 个参数，增长来自 encoder 输入层、denoiser 输入层和输出层；主体隐藏层不随 gene 数增加。权重本体不包含 optimizer、梯度、EMA 或激活显存。

18k 实际五折记录：每折约 56.5 分钟，五折共 16,965 秒（4.71 小时）；每折约 840 个 test wells 的 1,000-step DDIM 预测耗时约 42–44 秒。设备总显存采样峰值约 7.6 GB，包含约 2.5 GB 桌面基线，不能当作精确进程峰值。来源：18k `models/squidiff/README.md`、`results/run_cards.csv`。

因此“扩到 15k 必然需要大改架构、参数增加一个数量级、超过既有计算预算”的淘汰前提不成立：现有 18k 配置已经完成运行，参数增长约 3 倍，且在原 16 h/model 预算内。这不证明能够以相同成本修复收敛。

**后续决策（2026-10-26，用户已确认）：** 基于 5k/18k 两轮均未收敛、继续排障收益不明确的机会成本，将 Squidiff 排除出后续选型与训练/排障投入。保留历史代码、环境与全部实验结果，历史注册表不变；未来比较只包含其余六个候选。见[当前选型决策](model_selection_decisions.md)。下文对 Squidiff 的诊断建议仅在用户重新开启该候选时适用。

预算公平不等于 epoch 数相同，也不等于作者超参数原样搬来就适合新域。ADR 0027 明确限定 random init、禁止 HPO、单种子，并删减外部预训练。该合同回答的是“受限资源下的从零训练比较”，不是“各方法在合法资源与充分调优下的最佳表现”。这些约束有可复现性价值，但对依赖预训练的方法影响不对称。

外部纯化学结构预训练并不自动构成测试 CM response 泄漏。未来若允许，需要区分结构先验与药物响应先验，审计药物/上下文重叠，并在新协议中给出相同资源许可，而不是直接加载任何公开 checkpoint。

## 3. 数据独立性与泛化声明需收窄

**5k HVG 使用了 test wells 的表达分布。** `scripts/build_gene_set.py:1-18` 与 `REPORT.md:52-55` 明确使用三个 source 的全部 wells。它不使用显式 Δ 标签，但输入就是待预测的处理后表达；因此不能同时宣称“测试表达完全没有参与任何预处理”。这是已声明的 transductive feature selection，影响幅度未量化，不代表 PRnet 优势必然消失。更严格的确认实验应使用外部冻结基因面板，或仅在外层训练集内选基因。

**connectivity-block 不等于 scaffold 或 similarity holdout。** `scripts/build_data.py:521-523` 基于 InChIKey 前 14 字符建组，实际为 460 drugs / 458 blocks；能防相同连接结构及部分立体异构体跨角色，不能阻止相似骨架/近邻类似物跨角色。ADR 0027 的“structural analogues cannot leak”表述过宽。现有结果是 cold identity，不是 cold scaffold。

**共享板与对照限制外推范围。** 各 fold 仍覆盖同一 16 plates；共享同板 DMSO 对既定任务合理，但不是 cold plate、cold batch、cold cell context。每 drug-dose 只有约 3–6 个汇总 well，不能等同于作者大规模单细胞分布训练。

**重复观测测试结果影响最终选择。** 5k 是看过 18k 后追加的敏感性轴；后续修复与调参会继续消费现有测试证据。修复复算可以做，但最终确认不能继续把这组数据当完全未见的测试集，应引入新的独立验证或嵌套选择。

## 4. 指标一致，但含义与模型选择仍有限制

- `metrics.py` 实际使用 `average_precision_score`，即 AP；不是梯形积分 PR-AUC。“@200”表示真值 top/bottom-200，AP 排名范围仍为全部基因，不是只评价预测前 200。
- 真值由观测 Δ 排序定义，不包含显著性检验、FDR 或重复间一致性门槛；宜称方向 top/bottom-k 恢复，不直接称“显著 DEG 恢复”。弱响应组也分配 top-k，需要补充重复一致性/信号强弱分析。
- 5k 的正例比例为 4%，18k 为约 1.10%；绝对 AP 不直接横比。即使比较同轴相对 baseline 的 gap，也无法单独分离“训练维度改变”与“评价基因集/真值 top-k 改变”。
- PRnet 在 5k 的 AP 优势不等于幅值最好：报告 R² −0.2305、RMSE 0.1731，delta_mean 为 0.0164、0.1657。下游若使用 Δ 幅度，必须单独验收。
- 强 delta_mean 提示共性响应贡献较大；补充药物身份置换/消融、相对均值响应的残差指标、剂量分层、弱/强响应分层，检查模型是否学到足够的 drug-specific 信息。附加分析需标明 exploratory，不能事后更换主指标选赢家。
- checkpoints 的验证目标并未统一为部署任务：PRnet 为 control→treated NLL；chemCPA 验证中输入 treated 并重建 treated；CPA 原生 `cpa_metric = r2_mean + 0.5*r2_var + exp(disnt_after-disnt_basal)`（`pert_model/CPA/source/cpa/_task.py:599-629`，`_model.py:545-565`）。这不属于 test leakage，但与 control→cold-drug Δ 的最终目标存在错位。应明确保留作者选择准则还是增加独立 valid counterfactual 选择分支。
- 只有单训练 seed；drug bootstrap CI 是对固定已训练模型的测试药物重采样，不包含训练 seed、split、面板选择与多模型择优的不确定性，也没有化学相似性聚类修正。
- `ideal_fit=full` 是输入形式适配标签；`winner_eligible` 仅查五折 completed。它们均不检查收敛、配方完整、输出尺度正确；Squidiff 仍可被标 full/eligible。需要分开报告运行完成、语义正确、训练正常、任务适配。

## 5. 推进前的最小闭环

1. 先修复并验收 PerturbNet 5k 输出空间，复算受影响结果；冻结原始结果，不覆盖历史档案。
2. 对候选及关键竞争者增加语义验收：数据→模型→评测空间 round trip、全轴/子轴 library 守恒、drug/dose 条件响应、对照重建、tiny-subset overfit；Squidiff 先在其有效维度证明能训练/采样。
3. 新协议区分“统一资源从零训练”与“完整方法可用资源”问题；训练内部可有模型特定预处理，但外部评价目标、信息边界与输出转换必须一致。允许合理且事先限定的 validation 调优，保留预算/搜索次数/收敛证据。
4. 优先复核 PRnet、修复后的 PerturbNet、CPA、完整或明确受限的 chemCPA；保留 delta_mean。多 seed，并使用独立确认集或嵌套评估，避免继续在旧 test 上选最优配方。
5. 若目标仍是 whole-transcriptome MIRI，额外验收 5k 外的输出覆盖和下游基因集覆盖；目前 5k 领先不能自动满足 18,164 基因的原始目标。

建议当前决策记录为：**PRnet 暂列优先候选；5k PerturbNet 结果需要修复复核；尚未完成足以证明各方法充分发挥的最终模型选择。**
