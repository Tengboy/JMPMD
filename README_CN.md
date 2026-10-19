# JMPMD 正文实验开源复现包

本包整理当前博士论文正文使用的 **7 组图**：A1/A2 的角度、定位和检测三组图，A4 收敛图，以及 B1 导频数、A3 BS 天线数、B3 RIS 边长三组资源扫描图。QPSK 和 Gaussian 均采用严格矩，与当前理论口径一致。

本包服务于当前修订实验的复现。它包含相对原始 JMPMD 代码的实现补全、数值修正和初始化扩展，不能称为原论文所有结果不加修改的逐点复刻。详细说明见 [方法与模型](docs/METHODS.md)。

## 三种使用方式

### 1. 直接查看或重绘已经发表的结果

PDF/PNG 位于 `results/published/figures/`。主曲线数据保存在 `summary.csv`；A4 四条曲线的完整 200 轮统计在 `convergence.csv`，绘图时显示前 50 轮。

```matlab
cd('解压后本包所在目录');
setup_jmpmd;
plot_results;
```

默认导出至 `replotted/`，生成 PDF、300 dpi PNG 和可编辑 FIG。图内无图名，子图编号在下方，中文轴标和图例保留算法缩写；QPSK 理想 CSI 仅一条完整协方差联合 MAP 曲线。正文图注和 LaTeX 引用片段见 `docs/FIGURES.md`、`docs/figure_snippets.tex`。

绘图需要中文字体。Windows 优先宋体，Linux/macOS 可安装 Noto Serif CJK SC，或指定 `plot_results([],[], 'Font','你的中文字体名称')`。已导出的 PDF/PNG 不依赖使用者本机的 MATLAB 字体。

### 2. 验证安装与运行流程

```matlab
validate_release;       % 独立矩公式验证、60 个历史样本对照、恢复及失败保留检查
run_experiments;        % 默认为 smoke：全部220配置，每配置2次，最多4轮
verify_run(fullfile(pwd,'runs','smoke'));
```

MATLAB R2023b 为已测试版本。串行使用基础 MATLAB；并行需要 Parallel Computing Toolbox。默认不建立并行池、不安装定时任务。小样本流程检查不用于解释算法性能，QPSK 少量样本中的零 SER 也不会被人为替换成正数以绘制对数图。

### 3. 从头复算正文实验

```matlab
run_experiments('Profile','full','Workers',12);
verify_run(fullfile(pwd,'runs','full'));
plot_results(fullfile(pwd,'runs','full'),fullfile(pwd,'figures_full'));
```

完整运行包含 **220 个去重配置、130 万次估计器试验、232 个视图**：QPSK 75 万次，Gaussian 55 万次。A4 的 Q/G、0/6 dB 共 2 万条轨迹与对应 A1 条件共用样本。44 个条件需要参考量，合计 26 万条，全部从头计算，无须历史参考缓存。

完整运行较耗时；本次打包验证没有重新运行这 130 万次试验。默认 `smoke` 用于先确认安装无误。完整复算需要比发布包大得多的磁盘空间，建议输出位置预留数 GB；并行 worker 的内存也需单独考虑。

可分组运行，例如先复算 RIS 边长扫描：

```matlab
out=fullfile(pwd,'runs','B3');
run_experiments('Profile','full','Groups',{'B3'},'Workers',4,'Output',out);
verify_run(out);
plot_results(out,fullfile(pwd,'figures_B3'));
```

也可用 `'Kinds',{'Q'}` 或 `{'G'}` 只算一种调制。当前正文图为 Q/G 合图，缺另一种数据时绘图器会说明缺失，不会用随包结果混入新运行结果。

## 中断与恢复

每批结果先写临时文件，再提交为不可覆盖的 `chunks/c*.mat`。重复执行**相同参数与 Output** 会扫描已提交批次，仅运行未提交样本。修改源码、实验选择、运行档位或 MATLAB 版本后，恢复入口会拒绝混算；需要换新的输出目录。

主动分段运行可设置 `'MaxBatches',10`，或在输出目录放一个 `STOP_AFTER_CHECKPOINT` 文件；删掉停止标记后继续。若电脑关机，先确认 `process.json` 对应主进程已经退出，再人工处理遗留 `RUNNING.lock`；不能凭时间戳陈旧或 bootstrap PID 判断。详见 [运行协议](docs/REPRODUCING.md)。

## 实验与统计约定

- A1/A2：−5:1:10 dB；QPSK 7–10 dB 每配置 15000 次，其他条件均 5000 次。
- 默认 3 用户、3 RIS、BS 32 天线、RIS 10×10、帧长 40、导频 10、位置先验协方差 0.1I。
- A3：BS 16:4:40，SNR=6 dB，位置先验中心初始化。
- B1：导频 5/10/15/20，总帧长固定 40，SNR=6 dB。
- B3：RIS 边长 8:2:16，发射功率 25 dBm、噪声 −125 dBm 固定，使用局部导频角度初始化；其 SNR 随阵列变化。
- 全部最多 200 轮，连续两轮相对变化小于 1e−5 才判收敛；A4 画前 50 轮，提前停止的轨迹按末值延拓。
- 主指标保留失败、NaN/Inf、局部坏解与极端尾部；不得删除后再报告“改善”。均值、SEM 与调用/整次失败有独立列。
- W 是无定位反馈的迭代估计，结束后做一次角度定位；正文定位图不展示该补充定位曲线。
- A5/A6、RIS 数量 B2、严格/近似矩对照 B4、旧 Figure 3 和求峰器消融不属于本包正文范围。

## 数据和代码组织

`src/` 保留六个估计器的更新顺序与数值逻辑，`experiment_catalog.m` 定义实验矩阵，`simulate_trial.m` 生成共同观测并计算估计器/参考量。新控制器负责保存、恢复和统计。`metadata/source_provenance.json` 记录原始来源及 SHA-256，`tests/golden_samples.csv` 为少量回归样本。

**发布包不含原始实验 MAT、检查点、参考缓存、MATLAB 偏好或并行池缓存。** `results/published` 为只读基准；模拟运行后本地产生的 MAT 用于恢复，不应提交到开源仓库。随包 `SHA256SUMS.txt` 可检查文件完整性。

许可采用 [MIT](LICENSE)。原作者/贡献者的署名与来源说明保留在 [NOTICE.md](NOTICE.md)，引用说明见 [CITATION.md](CITATION.md)。
