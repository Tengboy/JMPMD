# Release validation, 2026-09-30

The archive was extracted to a separate directory containing a space in its path. A fresh MATLAB R2023b session used the default MATLAB path plus only the extracted release. Historical experiment/source directories were not placed on that path. The following checks passed:

| Check | Evidence |
|---|---|
| Numerical dependency closure | Base MATLAB and optional Parallel Computing Toolbox; no external project files required |
| Strict steering moments | Independent 1024-point circular quadrature on a rectangular 3×2 array |
| Gain/data Gram | Joint finite-state enumeration, complex off-diagonal covariance, full noise precision, Hermitian/PSD and deterministic-limit checks |
| Historical regression | 60 frozen full-iteration samples over 30 configurations, including all six estimators, NLoS, pilot counts, BS and RIS scans |
| Maximum scaled regression discrepancy | 4.2238610236926072e−15; time and memory excluded |
| Full A4 trace behavior | Four Q/G × 0/6 dB single trials, 200-row traces, correct final-position error and held final value |
| Full catalog smoke | 220 conditions, 232 views, 440 short estimator trials, 88 reference calculations, eight short traces |
| Smoke numerical status | Zero whole-trial, reference or primary nonfinite failures |
| Resume | One committed trial followed by continuation; original committed bytes unchanged; completed repeat adds no batches |
| Invalid resume | Changed experiment selection, existing lock and duplicate committed trial all rejected |
| Failure aggregation | Injected NaN trial retained in sample count, failure count and ordinary means |
| Published redraw | Seven main figures, 70 curves, 1090 points; each point equals the accepted thesis CSV/spec |
| Figure inspection | No internal figure titles; lower subplot letters, compact resource panels and padded Chinese legends; representative PNG exports visually checked |
| Distribution | No raw MAT, FIG, log, preferences or historical caches; text files contain no machine-specific drive paths |
| Frozen inputs | Original numerical source hashes remain unchanged |

The freshly generated smoke CSV was also plotted: five figures were exportable; two QPSK log-axis detection plots were explicitly skipped because two-trial ideal-CSI SER was zero. No pseudo-counts or published values were inserted.

`metadata/validation_results.json` contains the structured results and tested source signature. `metadata/published_data_audit.json` records exact equality to the selected thesis curves. `SHA256SUMS.txt` covers the distributed files; run `python tests/check_distribution.py` immediately after extraction to verify them, before generating new local outputs.

**Scope:** this is release/integration validation, not another full 1,300,000-trial Monte Carlo run. Published CSV/PDF/PNG originate from completed formal experiments. Linux/macOS and other MATLAB versions were not run for this release. Timing, memory and rare numerically sensitive local solutions can differ by platform.
