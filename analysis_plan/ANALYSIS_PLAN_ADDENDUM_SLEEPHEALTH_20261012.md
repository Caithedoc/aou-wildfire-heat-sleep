# 分析方案附录: 复合睡眠健康结局 (2026-10-12)

**状态**: 预设。本附录在任何复合结局模型运行之前定稿; 定稿后不回改定义,
新想法写进"事后探索"节。加入时点与理由如实写进稿件 Methods:
"added after coauthor review (Z.Q., 2026-10-09/10) and before any composite model was fitted"。

母方案: `ANALYSIS_PLAN.md` (2026-08-28, 主结局 TST, 预设次级 WASO / deep / REM / efficiency)
与 `ANALYSIS_PLAN_P2_MECHANISM_20260831.md` (时相 / 心率通路, 2026-10-12 并入 P1)。
本附录只新增**一个主复合结局与一个计分版**, 不改主结局, 不改任何已预设次级结局。

---

## 1. 为什么加

- 现有结果是分钟尺度的均值差 (复合夜 TST −4.03 min), Limitations 已承认没有 MID。
  复合结局把它翻译成"达标夜概率下降多少个百分点", 可换算成每千夜多出多少不达标夜。
- 时长、效率、连续性合成一个有共识定义的健康结局, 术语用 **sleep health**
  (Buysse 2014, doi:10.5665/sleep.3298), 不用 sleep quality (主观量表的词)。
- 预期两种结果都可发表 (见 §5), 本附录事先写明各自解读, 防止结果驱动叙事。

## 2. 定义 (v1, 第 1 步 VM 运行即用, 不需新数据)

每夜 (主睡眠段) 三个成分, 各记 0/1:

| 成分 | 变量 | 达标切点 | 依据 |
|---|---|---|---|
| 时长 | `minute_asleep` | 420–540 min (7–9 h) | Hirshkowitz 2015, doi:10.1016/j.sleh.2014.12.010 (成人 26–64 岁推荐 7–9 h) |
| 效率 | `sleep_efficiency` × 100 | ≥ 85% | Ohayon 2017, doi:10.1016/j.sleh.2016.11.006 (NSF 睡眠质量共识) |
| 连续性 | `minute_awake` (WASO) | ≤ 20 min | Ohayon 2017 |

- **主复合结局 `healthy_night`**: 三项全部达标 = 1, 否则 0。模型中按 0/100 编码,
  系数即百分点。
- **计分版 `sh_score`**: 三项之和, 0–3, 作补充。
- 任一成分缺失 → 两个复合结局均置缺失, 不以 0 代替。
- 切点统一用成人值, 不按年龄分层 (v1 的 parquet 无年龄列)。
- **预设敏感性 `healthy_night_s`** (附录用): 消费级腕戴设备相对 PSG 系统性高估 WASO、低估效率,
  因此把连续性与效率切点各放宽一档: WASO ≤ 30 min、效率 ≥ 80%, 时长不变。目的是检查 v1 基率偏低时
  结论是否依赖切点; 它不替换 v1, 也不因结果好坏而互换。

### v2 (时相重建 `04b` 之后, 在 `20_timing_joint.py` 中加入, 作敏感性而非替换 v1)
- 第 4 成分: 入睡潜伏期 `minute_to_fall_asleep` ≤ 30 min (NSF 公开摘要值; 全文将 ≤15 评为 good、16–30 为 appropriate, 开跑前按全文核对并在此处记录)。
- 第 5 成分: 觉醒次数 (`sleep_level` 中 ≥5 min 的 wake episode) ≤ 1 (≥65 岁 ≤ 2)。
- 年龄特异时长切点: ≥65 岁 7–8 h。
- v2 四成分 / 五成分复合与 v1 并列报告; v1 保持主复合结局。

## 3. 模型

与 `18_secondary_joint.py` 完全相同: person × 年月 × dow 分层固定效应, tmean 自然样条 5 结点,
四套规格 (four / int / marg / tlin), SE person / zip3 / zip3-day / 双向, 正文报 zip3-day, WCB 999 次。
- `healthy_night`: 线性概率模型 (与主结局同一估计量, 系数 = 风险差百分点); 不做条件 logistic,
  理由: 与主结局框架一致、可直接做 excess beyond additivity 的 delta-method 对比。
- `sh_score`: 线性, 单位 = 达标项数。
- 每千夜换算: 风险差 × 10; 复合夜归因不达标夜数 = 风险差 × both 夜数 (17,295) / 100。

## 4. 报告规则

1. 复合结局与三个成分**并排**报告 (同一表), 成分用各自原单位。这是 composite 的标准做法, 也是
   防止"由时长一项带着走"的唯一办法。
2. 只有 `healthy_night` 进正文主表; `sh_score` 与 v2 进附录。多重性: 复合结局计入现有
   "次级结局带 CI 不做校正" 的预设条款, 不单独加族。
3. 四格描述 (`rev7_secondary_cells.csv`) 必须报各格达标夜比例 (基率), 因为 Fitbit 的 WASO 偏高,
   `healthy_night` 基率可能偏低 (预期 20–40%), 要让读者看到分母。
4. 1–20 的计数一律不呈现 (AoU 披露规则)。

## 5. 事先写定的解读

- **复合夜 `healthy_night` 风险差显著为负, 且不止时长成分为负**: "sleep health 整体受损",
  标题可升为 "... and Sleep Health"。
- **风险差为负但完全由时长成分驱动 (效率、WASO 成分零)**: "损失集中在时长维度, 连续性未受影响",
  与已知的"少睡但没睡碎"一致, 支持时相/体温调节通路而非炎症-觉醒通路。标题保留 sleep duration。
- **风险差为零而 TST 均值差仍为 −4 min**: 均值位移未跨过达标切点, 说明效应是整体分布的小幅平移而非
  尾部恶化; 如实报告, 不换切点找显著。
- 三种情形都不改变主结局结论, 也都不触发对 v1 定义的回改。

## 6. 事后探索 (运行后才可填, 开跑前留空)

(空)

## 7. 脚本与产出

- v1 已写入 `workbench/18_secondary_joint.py` (OUTS 追加 `healthy_night`, `sh_score`, `healthy_night_s`; run_meta 追加基率)。
- 产出沿用 rev7 契约: `rev7_secondary_coefs.csv` / `_wcb.csv` / `_cells.csv` / `_bh.csv` / `rev7_run_meta.csv`。
- 锚点: `minute_asleep` four_M3 excess (zip3_day) 应复现 −1.78 (−3.13, −0.43)。
