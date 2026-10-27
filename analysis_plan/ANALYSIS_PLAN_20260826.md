# 野火烟 × All of Us 客观睡眠:执行方案 (2026-08-26)

工作题目: Machine learning-identified vulnerability to wildfire smoke-associated sleep disruption: a within-person analysis in the All of Us Research Program

前提已确认: AoU Researcher Workbench 已注册可用;Controlled Tier 提供真实未偏移日期 (Registered Tier 日期按人随机后移 1-365 天,不可用于逐日挂接)。

**Phase A 已完成 (2026-08-28, 结果见 RECON_RESULTS_20260828.md)**: CDR v9 (C2025Q4R6) 已挂入 workspace `aou-wildfire-sleep`;62,911 人 / 3,930 万主睡眠夜, 数据至 2024-12-02;2023-06 东岸事件在窗 (事件月 36,245 人在戴);zip3 覆盖 99.99%;中位佩戴 365 夜。主分析窗定 2016-01-01 至 2023-12-31 (smokePM v2.0 暴露端封顶)。
选刊路径 (journal-fit 2026-08-26): stretch = Nature Health (先 presubmission enquiry);Tier 1 = JAMA Network Open;Tier 2 = EHP / Environ Int;Tier 3 = Environ Res。

## 设计一页纸

- **Estimand**: 个体内急性效应。同一个人,烟日 vs 自身对照日的客观睡眠差,滞后 0 到 7 天恢复轨迹。
- **主设计**: time-stratified case-crossover (同人 × 同月 × 同星期几取对照日),辅以个体固定效应分布滞后模型 (DLNM, lag 0-7)。
- **暴露**: zip3-日 烟特异性 PM2.5 (Stanford ECHO Lab 产品,Childs 2022 EST, DOI 10.1021/acs.est.2c02934;v1 覆盖 2006-2020,v2.0 更新版覆盖至更近年份,下载时以实际末端为准)。主定义烟日 = 烟 PM2.5 > 20 µg/m3;敏感性 >10、>35。
- **结局** (逐夜,主睡眠段): 总睡眠时长、睡眠效率、WASO、入睡时间偏移、deep/REM 分钟数;**佩戴/未佩戴本身作为第二结局**。
- **关键共暴露**: 逐日温度必须同时进模型 (热与野火同季共现,热 → 睡眠效应已被 Liao 2025 证实,不调整温度 = 审稿必死)。温度用 PRISM 或 Daymet zip3 日值。
- **ML 三承重点** (逐相推进,先出干净的主效应再上 ML):
  1. 反事实预期睡眠基线 (个体历史 + 季节 + 星期 + 温度,梯度提升;锚 Heffernan 2025 AJE ML-CITS, DOI 10.1093/aje/kwae171)
  2. causal forest 脆弱性画像 (输入:年龄、慢病、剥夺指数、基线睡眠表型、空调可及性代理;GATES 校准检验防过拟合)
  3. 佩戴概率逆概率加权 (回应 Hurwitz 2026 佩戴偏倚框架, DOI 10.64898/2026.03.06.26347799)
- **样本参考量级** (来自已发表文献,以 recon 实测为准): Fitbit 睡眠 8k-23k 人 (BYOD 段) 至 59k 人 (设备发放段并入后, Patten 2026 Nat Med);烟事件受影响子集预计 2k-5k 人 (西部) + 2023-06 东岸事件若在数据窗内则样本大增。

## 执行相位

### Phase A: CDR recon (今天可跑, `workbench/01_cdr_recon.py`)
在 Workbench 里回答五个问题,全部有代码:
1. 工作区是 Registered 还是 Controlled Tier?CDR 版本号?(真实日期与 zip3 都要 Controlled)
2. sleep_daily_summary 的日期范围?**2023-06 东岸加拿大烟事件在不在窗内**(在 = 头号自然实验,东北部 AoU 参与者密集)
3. 有主睡眠记录的人数、人夜数
4. zip3 可得性与人数分布
5. 佩戴天数分布 (给 IPW 打底)

### Phase B: 队列与结局提取 (`workbench/02_sleep_cohort_extract.py`)
逐夜睡眠 + 人口学 + EHR 基线慢病 + SDOH 调查 + zip3,存 workspace 内中间表。

### Phase C: 暴露准备 (本地跑, `local_exposure_prep/03_smokepm_zip3_prep.py`)
本地下载 ECHO Lab zip 级烟 PM2.5 + PRISM/Daymet 温度 → 人口加权聚到 zip3-日 → 产出一张 csv 上传进 workspace。外部公开数据允许上传。

### Phase D: 挂接与主分析 (`workbench/04_casecrossover_analysis.py`)
zip3-日挂接 → case-crossover + DLNM → 主效应表 + 滞后曲线。

### Phase E: ML 三件套
D 的主效应稳了再上:反事实基线 → causal forest HTE → 佩戴 IPW 全线敏感性。

### Phase F: 负对照与稳健性
- 时间负对照:烟季前同月非烟年
- 结局负对照:与烟无机制关联的指标
- 设备型号/年份分层 (算法漂移)
- 热浪日剔除敏感性
- 细则见下面两节:数据清理规范 + 敏感性矩阵

## 数据清理规范 (2026-08-27 增补, Phase B/C 落地; 本节为预设, 开跑后不回改)

### 夜级过滤 (Phase B)
1. 只保留主睡眠段 (is_main_sleep 或 recon 实测等价字段)。小睡另存一表备用, 不入主分析。
2. 生理可信窗: 总睡眠时长 <3 h 或 >14 h 的夜剔除; 睡眠效率 >100% 或任何负值分钟数视为设备错误剔除。阈值是约定不是金标准, 用放宽窗 (1-16 h) 做敏感性。
3. 同人同日多条主睡眠记录: 保留 minute_asleep 最长一条; 完全重复行去重。
4. 夜的日期归属预设写死: 以睡眠起始日挂接暴露日。起始日 vs 醒来日错位列入敏感性。

### 人级纳入 (Phase B)
1. 纳入: 年龄 >=18, 有 Controlled Tier zip3, 窗内 >=60 个有效夜。夜数门槛做 >=30 / >=180 梯度敏感性。
2. 窗内 zip3 变更者 (搬家) 主分析剔除; 敏感性加回并用当期 zip3。
3. 睡眠障碍与用药人群不剔除, 改分层敏感性: EHR 里 OSA 诊断、失眠诊断、催眠药处方各一层。夜班工人若 SDOH/职业字段可识别, 做剔除敏感性。

### 佩戴与缺失 (Phase B, 喂 Phase E)
1. 佩戴日定义按 recon 实测字段定 (常见: 当日有任意心率或步数记录)。
2. 佩戴缺失不当 MCAR 处理: 烟日可能改变佩戴行为, 所以佩戴/未佩戴本身是第二结局 (方案既有), 同时给 IPW 建模 (Hurwitz 2026, DOI 10.64898/2026.03.06.26347799)。
3. 主分析 = complete nights; IPW 加权版为平行呈报的主敏感性, 两者同表并列。

### 暴露侧清理 (Phase C, 本地)
1. smokePM 的 0 与缺失必须区分: 产品对无烟日赋 0, 先核对覆盖掩膜与版本末端 (Childs 2022, DOI 10.1021/acs.est.2c02934), 再写缺失规则。
2. zip5 到 zip3 人口加权: 用 2020 Census ZCTA 人口, 权重年份写死; 2010 权重做敏感性。
3. 温度 (PRISM 或 Daymet) 同法聚 zip3-日; 两个温度源互为敏感性。
4. 挂接后做一致性核查: 抽 zip3 若干, 烟日序列对照已知 2018/2020/2023 事件日历, 对不上即停。

### Attrition 表
每条过滤规则一行, 报人数与人夜数瀑布。导出遵守 <20 禁令, 小格并入大类。

## 敏感性分析矩阵 (Phase D/F 执行; 每族回答一个具体威胁)

| 族 | 威胁 | 分析 |
|---|---|---|
| 暴露定义 | 烟日阈值任意性 | 主 >20 µg/m3; 敏感性 >10、>35; DLNM 连续暴露反应曲线为无阈值版本 |
| 暴露错分 | zip3 聚合与搬家 | 剔除 vs 纳入搬家者; PRISM vs Daymet 互换; 限高覆盖 zip3 子集 |
| 混杂 (共暴露) | 热与烟同季共现 (Liao 2025, PMID 41289815) | 温度不调 vs 线性 vs 非线性样条三档; 剔除热浪日; 热 x 烟交互项 |
| 设计规格 | 对照日方案与滞后结构 | time-stratified 主案 vs 对称双向对照; lag 0-7 主案 vs 0-14; DLNM 节点数变动 |
| 选择与测量 | 佩戴偏倚、算法漂移、志愿者结构 | IPW vs complete nights; 设备型号/年份分层; 夜数门槛梯度; 佩戴作为结局的缺失机制分析 |
| 证伪检验 | 设计本身失效 | 未来暴露负对照 (lag -1/-2 的烟不应预测今晚睡眠); 时间负对照 (非烟年同月); 结局负对照一项 (候选待与烟机制脱钩论证后定, 不硬凑) |

补充三条:
1. 多重性: 主结局预设为总睡眠时长, 主暴露窗预设为 lag 0-1 累积; 其余结局与滞后一律标 secondary, 全部带 CI 呈报, 不做正式校正 (环境流行病学惯例, 文中明说)。
2. 推断: case-crossover 用条件 logistic / 个体固定效应, 按人聚类稳健 SE。
3. 定量偏倚: 主效应报 E-value; 若佩戴 IPW 与 complete nights 结论相左, 升级做简易 QBA 而不是挑好看的报。

本节连同全文 = 预设分析方案。开跑前存档一份带时间戳 (OSF 或 workspace 内 commit), 当轻量预注册。

## 纪律红线

1. **披露规则**: 任何导出的表格计数 <20 一律禁止;分析全程留在 workspace,只导出聚合估计。
2. **代码里所有表名列名先跑 discovery cell 再用**,AoU 数据字典随 CDR 版本变,方案里的字段名一律以 INFORMATION_SCHEMA 实测为准。
3. **不做患病率或全国代表性表述** (非概率样本,无权重)。
4. 结果出来之前不联系期刊;效应平淡则 Nature Health 从路径删除,从 JNO 起投。
5. 零 em-dash;引用 ID 全部现场核实过,新增引用照此办理。

## 与综述线的分工

综述 (已立项,pitch 待发) §2.2 点名"野火烟 × 客观睡眠零大队列证据" = 给本篇插旗;本篇 Discussion 回引综述。两线并行,互不阻塞。

## 开工前 Cai 需要拍板的两件事

1. Workbench 里用 Python 还是 R?(脚手架先按 Python + pandas_gbq 写,R 版可再出)
2. 作者团队:建议 Cai 1st/corr + Geiger senior (Team G 变体),烟暴露方法如需可加一位空气污染建模合作者
