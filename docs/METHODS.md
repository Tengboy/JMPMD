# 方法、实现变化与统计口径 / Methods and retained assumptions

本包使用当前正文实验的严格矩版本。严格矩指在**既有均值场因子分解下完整计算所需二阶期望**，不意味着求得真实联合后验、全局最优解或证明全局收敛。以下变化分为实现补全、数值修正、初始化扩展和可复现工程，不把原论文已包含的期望写成新提出的理论。

## 1. 系数更新的二阶矩

令第 m 个 RIS、第 t 个时隙、第 k 个用户的投影响应为 ψ_mtk=ω_mt^T a(θ_mk)。角度因子按用户/轴独立，记其均值 p_mtk=E[ψ_mtk]、功率 d_mtk=E[|ψ_mtk|²]。符号向量 x_t 的后验均值为 μ_t，复协方差为 C_t=E[(x_t−μ_t)(x_t−μ_t)^H]。

系数更新需要 J_m=Σ_t E[A_mt^H A_mt]，其元素为：

\[
[J_m]_{ij}=\sum_t (\mu_{ti}^{*}\mu_{tj}+[C_t]_{ji})
\begin{cases}
d_{mti},&i=j,\\
p_{mti}^{*}p_{mtj},&i\ne j.
\end{cases}
\]

不同用户角度后验独立，因此非对角的角度期望可以分解；符号的用户间后验协方差不因此消失。复协方差的下标为 ji，代码用非共轭转置 `.'`。把软符号 μ 当确定值会漏掉这些协方差项。

原始代码已有角度响应功率矩的计算，不能概括为“原算法完全未考虑角度不确定性”。Gaussian 的具体旧入口还存在角度对角替换未启用的分支，当前统一使用 d。QPSK 是单位模，E|x_k|²=1，其旧系数 Gram 对角直接采用 d 时已经包含正确常模二阶矩；本次 Q 的关键补全包括非对角符号互协方差，而不是宣称常模对角也必然漏方差。

已知导频 C_t=0；确定角度 d=|p|²，退化为普通设计矩阵 Gram。源码：`src/functions/vb_gain_gram.m`、`vb_steering_moments.m`、`steering_power.m`。

## 2. 数据更新的完整信道二阶 Gram

每个 RIS 的系数后验向量保留完整用户协方差 C_ρ,m。令 P 为投影噪声的完整精度矩阵，H̄_t=E[H_t]，h̄_mt 为其第 m 行。独立 RIS 因子下：

\[
E[H_t^H P H_t]=\bar H_t^H P\bar H_t+
\sum_m P_{mm}\left(E[h_{mt}^Hh_{mt}]-\bar h_{mt}^H\bar h_{mt}\right).
\]

同一 RIS 行的二阶矩为 `(conj(μρ)*μρ.' + Cρ.')` 与角度二阶矩的逐元素积。修正包含系数非对角协方差，P 的跨 RIS 项也保留；没有把投影噪声擅自对角化。Gaussian 数据后验使用该二次型和单位复高斯先验。QPSK 枚举 4^K 个联合状态形成后验及完整符号协方差；SER 用联合 MAP，软 MSE 用后验均值。

原 Gaussian 数据更新已包含部分对角功率修正，本次是统一完整二阶量，而非所有更新都从一阶量开始。Hermitian 对称化只清理舍入误差，不剪负特征值。源码：`vb_data_gram.m`、`discrete_posterior_moments.m` 和六个 `estimate_*`。

## 3. 首次系数初始化与角度局部更新

JMPMD 和无定位反馈方案的首次系数估计使用已知导频，降低零/较差数据软估计导致的局部坏解；P 本来只在导频上估计。仅增加最大迭代轮数不保证消除局部坏解。

`Heuristic3.m` 按当前角度均值做局部上升，限制最大相位步长并回溯接受不降低目标的步长。旧的全周期找峰回退可能跳向远处旁瓣，当前移除这种回退。von Mises 曲率近似、位置反馈和其他原模型定义保留。源码不是全局求峰器，也没有保证所有初始化均落入正确主峰。

收敛判据只使用估计状态的变化，不使用真值误差：对状态字段取当前与上一轮差值范数，除以两者范数与 realmin 的最大值；所有字段的最大相对变化连续两轮小于 1e−5 才停止，上限 200。不同方法使用其适用状态字段。

## 4. B3 的先验局部导频角度初始化

RIS 边长变大，主峰变窄，位置先验均值对应角度的偏差可能导致坏局部解。B3 在进入变分循环前，用导频观测做一次先验约束的角度搜索：位置先验经局部线性化给出角度中心与方差；搜索限制在约 ±3σ 内，先粗后细三层、三轮坐标扫描；消除其他用户导频列的子空间，并剖面化未知复增益。

接口只传 `X_p, Y, A_BS, Omega, effective_noise_power, User_pos_c, RIS_pos, RIS_dir_x, RIS_dir_y`。不传真实用户位置、真实角度、真实系数或未知数据。只替换初始角度均值，角度先验和初始集中参数保持原定义。所有 B3 的 Q/G、J/P/W 都使用该初始化。

该搜索是 B3 的额外初始化步骤，不能把性能收益全部归因于严格矩。A3 BS 扫描仍用位置先验中心初始化。源码：`src/pilot_initial_angles.m` 及各估计器的可选 `initial_theta` 输入。

## 5. Baselines

- **JMPMD (J)**：角度、系数、符号与位置信息迭代更新。
- **仅导频方案 / pilot-only (P)**：参数估计只使用导频，完成后检测数据；不能把其变化解释为引入未知符号协方差参与导频参数估计。
- **无定位反馈方案 / no position feedback (W)**：估计过程中没有位置定位反馈。结束后从角度做一次最终定位；该定位结果属于补充分析，正文定位图未使用 W 曲线。
- **理想 CSI / ideal CSI**：QPSK 已知真实 H，保留完整投影噪声 C=σ²(A_BS^H A_BS)^{-1}，枚举联合符号并最小化残差 C^{-1} 二次型。等概星座下是联合 MAP/ML。它优化联合状态判决，不是普适逐符号 Bayes SER 下界。历史对角 C 度量只在 CSV 中保留，不绘为第二条基准。

## 6. 观测模型、SNR 与 NLoS

位置先验均值为场景中给定的三用户中心，真实位置按 0.1I₃ 高斯扰动。RIS 探测相位随机。默认 T=40、Tp=10；与原稿部分描述的 T=50 有区别。QPSK 单位模，Gaussian 为单位功率圆对称复高斯。

A1/A2 的 SNR=10log10(mean|Y_LoS|²/σ²)，噪声方差由 LoS 功率设置。加入 NLoS 不改变该噪声或 SNR 分子。B3 固定 Tx=25 dBm 和 noise=−125 dBm，其 ActualSNR 随阵列和样本变化，不是全程 6 dB。

NLoS `mixed20` 是当前实验定义：两跳级联包含 LoS×NLoS、NLoS×LoS 和 NLoS×NLoS，另有绕过 RIS 的单跳散射支路。级联与单跳模板混合后，按**BS 接收端投影前总干扰功率**统一缩放至 κ=10/15 dB；合成时保留交叉项。投影后的 κ 不保证等于输入 κ。估计器知道总干扰功率并将其加入有效噪声方差，但不建模其空间结构；这是实验的已知干扰功率假设，不能泛化为完全未知干扰下的结论。

## 7. References and error metrics

BCRB 标签沿用作者确认的项目定义：每次生成真值下的信息矩阵、加入指定先验精度，求逆后对各次界平均。**E[J^{-1}] 不一般等于 (E[J])^{-1}。** QPSK 定位/角度参考假设全部符号已知；Gaussian 为包括未知符号及其先验的联合信息条件。LoS 界不作为 NLoS 曲线的通用下界，也不跨调制强行作相同信息条件比较。稀疏 FIM/Kronecker/延迟结构加速保留原信息定义和坐标顺序。

- Theta：两个空间频率轴 `π × direction cosine` 的平均平方误差，单位 rad²，不是物理角度的度数 MSE。
- Position：每用户三维欧氏距离平方，跨用户平均，单位 m²。
- QPSK Data：联合 MAP 符号错误率 SER，不是 BER。
- Gaussian Data：未知符号后验均值 MSE。
- A4 X_MSE：Q/G 均为后验软符号均值 MSE；不是将 QPSK 的 SER 改名为 MSE。

## 8. Numerical fixes, diagnostics and release engineering

沿用已经复验的坐标/转置顺序、位置先验精度、方向响应延迟索引与 Taylor 定位导数修正；Taylor 协方差在返回位置计算。稀疏信息矩阵、结构化阵列矩和避免巨大中间矩阵属于实现加速，并未改变信号模型或界的求逆顺序。

整次异常、主指标非有限、到达迭代上限、未收敛、Taylor 中间调用失败和最终每用户超过 1 m 分别记录。Taylor 调用失败不等价于整次试验失败。A4 停止后按末值延拓，以同样的样本数报告每轮平均；活动样本数单独记录。均值不删除极端值或 NaN。SEM 是样本标准差/√N，不是重尾情况下的无条件可靠置信区间；未做多重比较校正。

发布包新增的部分是独立配置目录、串行/并行入口、可恢复检查点、从头参考计算和 CSV 重绘接口。六个估计器的数值更新顺序保留。去掉了 Q-P 中从未使用的 `gX` 对象构造和旧求峰器选择分支，减少外部代码依赖；60 个冻结样本用于检查这些封装没有改变正文算法结果。

## English summary

The release completes the required expectations under the existing mean-field factorization: symbol covariance in the gain Gram; full gain covariance and steering moments in the data Gram; full projected-noise precision. It retains the current local angular solver, pilot-based first gain update, and B3-only prior-local pilot angle search. It preserves the project's noise/NLoS definitions, information-bound ordering, initialization differences between A3 and B3, failure statistics and stopping rule. These are implementation completions and experiment extensions, not a new factorization or a global convergence guarantee.
