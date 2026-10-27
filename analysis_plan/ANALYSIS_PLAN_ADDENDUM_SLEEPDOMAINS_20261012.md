# 分析方案附录 2: 四域睡眠结局 (2026-10-12)

**状态**: 预设。Cai 2026-10-12 拍板 "先把睡眠做完", 方向 = "Liao 2025 镜像 + 时相"
(依据 `AOU_OUTCOME_LANDSCAPE_20261012.md` §G)。本附录在任何时相 / 结构域模型运行之前定稿,
定稿后不回改定义; 新想法写进 §7 "事后探索"。稿件 Methods 如实写明加入时点:
"secondary sleep domains were prespecified on 2026-10-12 after coauthor review and before any
timing or architecture model was fitted".

母方案: `ANALYSIS_PLAN.md` (主结局 TST; 预设次级 WASO / deep / REM / efficiency),
`ANALYSIS_PLAN_P2_MECHANISM_20260831.md` (H2 时相预测 B1–B4),
`ANALYSIS_PLAN_ADDENDUM_SLEEPHEALTH_20261012.md` (复合结局 healthy_night / sh_score)。
本附录不改主结局, 不改已预设的任何定义, 只把次级结局组织成四个域并补齐时相域。

---

## 1. 结局清单 (逐夜, Fitbit main sleep, 单位写死)

| 域 | 变量 | 定义 | 来源列 | 先例 |
|---|---|---|---|---|
| 时长 | `minute_asleep` | 主睡眠段睡着分钟 | sleep_daily_summary | Liao 2025 主结局 (**本稿主结局, 不动**) |
| 连续性 | `eff_pp` | 100 × minute_asleep / minute_in_bed (与 Fitbit sleep_efficiency 一致时直接用后者) | 同上 | Liao 2025 |
| 连续性 | `minute_awake` | 入睡后清醒分钟 (WASO); classic 与 stages 两种模式均为准确值 (AoU 文档) | 同上 | Liao 2025; Master 2025 |
| 连续性 | `n_awake5` | 主睡眠段内 ≥5 min 的 wake / awake 片段数 | `sleep_daily_summary_counts` 若存在, 否则由 sleep_level 计数 | P2 预测 B4; Master 2025 |
| 时相 | `onset_min` | 入睡时刻相对 sleep_date 当日 00:00 的分钟数 (前一晚 22:30 = −90) | `sleep_daily_summary_ext` 起始时间, 否则 sleep_level 首段起始 | Liao 2025 (相对当地午夜) |
| 时相 | `offset_min` | 末次醒来时刻相对同一 00:00 的分钟数 | 同上 (结束时间) | Master 2025 sleep offset |
| 时相 | `midsleep_min` | (onset_min + offset_min) / 2 | 导出 | Coleman 2026 |
| 时相 | `midsleep_shift` | 本夜 midsleep_min − 前一夜 midsleep_min, 仅当前一夜为连续日历日 | 导出 | Coleman 2026 (夜级规律性) |
| 时相 | `minute_to_fall_asleep` | 入睡潜伏期 | sleep_daily_summary | P2 预测 B1 |
| 结构 | `minute_light` / `minute_deep` / `minute_rem` | 各分期分钟 | sleep_daily_summary | Liao 2025; Elemosho 2026 (分钟) |
| 结构 | `pct_light` / `pct_deep` / `pct_rem` | 100 × 分期分钟 / minute_asleep | 导出 | Zheng 2024; Kang 2026 (占比, 敏感性) |
| 复合 | `healthy_night`, `sh_score`, `healthy_night_s` | 见 SLEEPHEALTH 附录, 不变 | 导出 | 无 AoU 先例 |
| 复合 v2 | `healthy_night_v2` | v1 三成分 + 潜伏期 ≤30 min + n_awake5 ≤1 全达标 | 导出 | SLEEPHEALTH 附录 v2 |

时区: Fitbit API 的起止时间为用户本地时间 (无时区标记)。开跑前先看 onset_min 的分布 (19 脚本 discovery
输出): 众数应落在 −150 至 +60 min (21:30–01:00); 若整体平移 5–8 小时则为 UTC, 届时按 zip3 → 时区
换算后再入模型, 并在本附录 §7 记录。

## 2. 夜集 (两套 N, 照 Liao 2025 分开报)

- **时长集**: 全部分析夜 (现 1,633,414 夜 / 30,268 人)。用于 minute_asleep、eff_pp、minute_awake、复合 v1。
- **结构集**: 分期分钟非空的夜 (= Fitbit stages 模式且有心率的夜; classic 模式夜分期为空)。用于
  三个分期分钟、三个占比、n_awake5、时相四变量、潜伏期、复合 v2。
- 两套都沿用主分析的夜级窗 (180–840 min, 效率 ≤100 %) 与层定义 (person × 年月 × 星期几, 层内 ≥2 夜)。
- 设备模式旗标 `mode_stages` = minute_wake 非空。进 cells 描述表, 并作分层敏感性 (与 src_id 分层并列)。

## 3. 模型

与 `18_secondary_joint.py` 完全相同 (正文 M3 规格): 层固定效应 + tmean 自然样条 5 结点 (全样本结点,
结构集沿用), 四套规格 four / int / marg / tlin, SE person / zip3 / zip3-day / 双向, 正文报 zip3-day,
int 规格的 sxh 做 zip3 层级 WCB 999 次。全部结局线性; 时相变量单位 = 分钟, 占比单位 = 百分点,
复合单位 = 百分点 (0/100 编码)。

## 4. 报告规则

1. 正文一张表 "Sleep domains": 行 = 时长 1 + 连续性 3 + 时相 4 + 结构 3 + 复合 1 (healthy_night),
   列 = 热单独 / 烟单独 / 复合 / excess (four_M3, zip3-day CI); 占比版、潜伏期、复合 v2 与敏感性入附录。
   表头注明每行用的夜集与 N。
2. 与 Liao 2025 Table 2 逐行对照放 Discussion 一段 (他们每 +10 °C 异常: TST −2.19 / −2.63 min,
   SE −0.03 %, WASO +0.05 min, onset +1.66 min); 本稿热单独格的方向与量级和它对表。
3. 多重性: 全部计入既有 "次级结局带 CI 不做正式校正" 条款; BH q 值按域分族仅作附录参考。
4. 1–20 的计数一律不呈现; 时相变量的格均值若由 <20 夜构成同样置空。
5. 术语: duration / continuity / timing / architecture / sleep health; 不用 sleep quality。

## 5. 事先写定的解读 (与 P2 方案 §3 对齐)

- 复合夜 onset_min 推迟且 offset_min 不提前或提前幅度更小 → "窗口由前端截短", 支持 H2 时相位移
  (体温散热 / 行为通路), 与 minute_awake 零结果一致。
- 复合夜 onset 不动而 minute_asleep 仍降 → 截短发生在后端或中段; 看 offset_min 与 n_awake5。
- n_awake5 与 minute_awake 同为零 → 碎片化通路彻底排除 (两个独立证据)。
- 分期分钟降但占比不变 → 整体缩短; 占比变 (REM 占比降) → 结构改变, 单独讨论。
- 任何结果都不触发 §1 定义回改。

## 6. 脚本与产出

- `workbench/19_sleep_domains_rebuild.py`: discovery (`_ext` / `_counts` / sleep_level 列名, onset 分布) →
  重建 `analysis_nights` (补全列 + 时相 + 觉醒次数) → 刷新 `analysis_nights.parquet`。
- `workbench/20_sleep_domains_joint.py`: 四域结局 joint 模型, 产出 `rev8_domains_coefs.csv` /
  `_wcb.csv` / `_cells.csv` / `_bh.csv` / `rev8_run_meta.csv`; 锚点 minute_asleep four_M3 excess (zip3_day)
  应复现 −1.78 (−3.13, −0.43)。
- 18 脚本保留不动 (复合 v1 已在其中); 20 脚本覆盖 18 的全部结局, 跑 20 即可。

## 7. 事后探索 (运行后才可填, 开跑前留空)

(空)

(运行后记录, 2026-10-12 VM rev8, 不改 §1–§6 任何定义)
- 数据字典实测: `sleep_daily_summary_ext` 与 `sleep_daily_summary_counts` 均存在, 时相直接取 sleep_start_time / sleep_end_time (本地时间, 入睡众数 22–23 点), 不走 sleep_level 重建。
- 复合结局基率远低于 §预期: healthy_night 1.33 %, healthy_night_s 2.44 %, v2 0.046 % (Fitbit WASO 均值 55 min)。按 SLEEPHEALTH 附录 §5 第三种情形处理: 如实报告, 不回改切点; 复合层降入附录。
- 入睡潜伏期 (minute_to_fall_asleep) 在 Fitbit 中均值 0.18 min, 无信息量, 列为方法学备注。
- 事后探索候选 (未运行, 仅记录): 以 WASO ≤ 第 25 百分位 (样本内) 替代固定 20 min 的复合定义; 时相域用 sleep_type 区分 classic / stages 分层。
