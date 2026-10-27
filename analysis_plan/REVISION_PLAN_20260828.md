# Compound-first 修订计划 (2026-08-28, 按 Cai 优先级)

依据 Cai 审阅意见 (7 项 + 1 个模型规范问题)。分析数据在 AoU Workbench VM,
本地只有聚合结果 → 所有新数字由 `workbench/11_revision.py` 在 VM 内一次跑齐,
产出 `rev_*.csv` 后本地完成图表重构。

## 状态总览

| # | 事项 | 状态 |
|---|---|---|
| P7 | Table 3 "Hot only 194,276" 数字问题 | **已处置**: 194,276/1,633,414=11.89% ≈ 正文 total hot 11.9%, 极可能是总热夜数误放 hot-only (真 hot-only ≈ 176,108); smoky-only 179,491 也对不上 12.09% (差 ~150 夜)。两格已改 "--" + Word 批注, `rev_cell_counts.csv` 精确回填 |
| 模型 | Methods "all models linear+quadratic temperature" 与 compound 模型不符 | **已改写**: 经查 06_phaseE 脚本, 交互模型只调线性 tmean, 四格模型无连续温度调整。Methods 已按实改写 + 批注 (VM 版本为准, 投前复核) |
| P4a | Table 2 只报 2/4 secondary outcomes (selective reporting 风险) | **半修复**: 04 脚本实际跑了全部 4 个 (WASO/deep/REM/efficiency), REM+efficiency 系数已在 VM 的 phaseD_results.csv, 只是没转录。取回即补全, 无需重跑 |
| P4b | Table 3 加 additive expectation / excess 行 | **已完成** (−2.29 / −1.73 已入表, 交互行原有) |
| P1 | continuous joint figure 取代四柱 Figure 2 | 待 VM: spline joint model + P50/75/90/95 曲线 → `rev_joint_curves.csv` |
| P2 | compound 专属 robustness (3×3 matrix + excl2023/东西/设备/IPW/宽窗) | 待 VM → `rev_matrix3x3.csv`, `rev_compound_sens.csv` (IPW 与宽窗需复用 06 E3 权重 / 09 重抽样本, 脚本内已标) |
| P3 | 新 Figure 1: flow + 年月暴露频次 + 区域分布 三 panel | 待 VM → `rev_fig1_timeline.csv`, `rev_fig1_region.csv`; flow 数字已有 (table_attrition) |
| P4c | Table 1 加 n(%)、race/ethnicity、region、device、ADI quartile、education、informative strata | 待 VM: `rev_table1_addons.csv` + BigQuery person/survey 拉取 (08_tables.py SQL 模式, n<20 抑制) |
| P6 | Table 4 / subgroup forest 加 formal P_interaction | 待 VM → `rev_subgroup_pint.csv` (需先导出人级 flag 表, 脚本内已标) |
| P5 | causal forest → Supplement; distributed lag / robustness forest → Supplement 或缩 panel | 待 P1-P3 数据到位后一并重排主图 (避免中间状态), 见下"目标主图结构" |

## 目标主图结构 (P1-P3 完成后)

- **Figure 1** = 三 panel overview: A flow (S1→S6), B 2016–2023 年月 smoke/heat/compound 夜数 (堆叠条/热图, 直面 2023 占 57%), C 区域分布
- **Figure 2** = continuous joint response curves (横轴 smoke PM2.5 0–60, 四条曲线 = 温度 P50/75/90/95, 带 CI; 或 contour surface 作 alt)
- **Figure 3** = interaction forest: subgroup | 点估计+CI | β (95% CI) / exposed nights / P_interaction 三列右栏
- → Supplement: 现 Fig 4 causal forest (calibration failed), 现 Fig 1 distributed lag, 现 Fig 5 sensitivity forest (或缩为新 Fig 2 的 inset), 3×3 matrix 表
- Table 3 重心 = observed / additive / excess / interaction (已改); Table 4 加 P_int 列

## 执行顺序

1. VM 重启实例 → 跑 `workbench/11_revision.py` (先补 person flag 表, 脚本有提示)
2. 取回 rev_*.csv + phaseD_results.csv (REM/efficiency 两行)
3. 本地: 重画 Fig 1/2/3 (600dpi PNG+SVG+PDF, 新编号), 重排主图/Supplement, 补 Table 1/2/4, 回填 Table 3 Nights
4. Methods 补 continuous joint model 描述 (P1 兼作温度规范敏感性), 更新 Word Count
5. 复核两条 Word 批注涉及的数字后删除批注

## 纪律

- 所有新导出遵守 n<20 抑制; 只出聚合系数/计数
- 每个新数字进稿前与 rev_*.csv 原文件比对, 不允许"顺手"手打
- 图表重排后 figure 引用顺序 = 编号顺序再核一遍 (上轮已理顺, 重排会再动)

---

## 执行完成 (2026-08-28 晚)

全部优先级已落实 (P1/P2/P3/P4/P5/P6/P7 + 模型规范), 细节见 submission_EnvironInt/BUILD_NOTES.md 末节。
唯一未跑: 宽窗 (60-960) 版 compound (需 CDR 重抽, 无缓存, 判定成本不值; 边际宽窗敏感性仍在 Table S3)。
VM 实例已停止。rev_*.csv 全部下载至 results_import/。
