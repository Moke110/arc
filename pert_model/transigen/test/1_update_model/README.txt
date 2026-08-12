# 1_update_model 测试目录说明

本目录用于验证升级后的 TranSiGen 核心模块（dataset/model/vae）是否可以完整执行一次
“读取数据 -> 模型预测 -> 输出结果”的流程。

## 目录结构
- input/: 测试输入数据
- output/: 模型输出结果
- run_test.py: 端到端测试脚本

## 输入文件（input）
- processed_data.h5
  - 作用: 测试表达数据主文件，包含 x1、x2、sig。
  - 说明: x1 为基线表达，x2 为扰动后表达，sig 为样本唯一标识。
- dataset_index.csv
  - 作用: 数据索引与映射文件，驱动 Dataset 逐样本读取。
  - 字段: LINCS_index, mid, cid。
- mid2smi.pickle
  - 作用: 分子整数 ID（mid）到 SMILES 字符串的映射字典。
- example_KPGT.pickle
  - 作用: SMILES到分子向量（2304维）的映射字典。

## 输出文件（output）
- predicted_profiles.h5
  - 作用: 保存模型预测的完整结果。
  - 内容: x2_pred（预测表达谱）、mid、cid。
- prediction_summary.csv
  - 作用: 预测结果的样本级摘要，便于快速检查。
  - 字段: mid, cid, pred_mean, pred_std。

## 运行方式
- 在项目根目录执行:
  - .venv/bin/python test/1_update_model/run_test.py
