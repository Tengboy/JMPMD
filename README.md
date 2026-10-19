# JMPMD: reproducible main-figure experiments

[中文说明](README_CN.md) · [Experiment protocol](docs/REPRODUCING.md) · [Methods and definitions](docs/METHODS.md)

This MATLAB release reproduces the **current thesis experiments** for joint multiuser positioning and data detection in a RIS-assisted system. It includes strict-moment QPSK and Gaussian estimators, their baselines, channel simulation, reference calculations, deterministic Monte Carlo runners, and the seven selected main-figure groups.

It is a revised research implementation, not a claim of bit-for-bit reproduction of every experiment in the original paper. The changes and retained assumptions are listed in [METHODS.md](docs/METHODS.md).

## Start here

Open MATLAB in the extracted release directory. MATLAB R2023b is the tested version. Base MATLAB is sufficient for serial execution and plotting; Parallel Computing Toolbox is optional.

```matlab
setup_jmpmd;
plot_results;                     % redraw the supplied published CSV; no simulation
validate_release;                 % math, 60 frozen samples, resume/failure checks
run_experiments;                  % smoke test: 2 trials/config, 4 iterations
verify_run(fullfile(pwd,'runs','smoke'));
```

The PDF/PNG files in `results/published/figures` are ready to view without MATLAB. Replotting uses Chinese labels, as in the thesis. It needs a Chinese-capable font; use `plot_results([],[], 'Font','Noto Serif CJK SC')` if needed. Replotting creates PDF/PNG/FIG under `replotted`, leaving the published results unchanged.

**Full simulation is a long computation, never the default.** It runs 1,300,000 estimator trials over 220 unique configurations and calculates 260,000 realization-wise references from scratch. The main plots use 232 views of those configurations; shared settings are executed once.

```matlab
run_experiments('Profile','full','Workers',12);  % optional parallel workers
verify_run(fullfile(pwd,'runs','full'));
plot_results(fullfile(pwd,'runs','full'),fullfile(pwd,'figures_full'));
```

For a smaller independent experiment:

```matlab
out = fullfile(pwd,'runs','B3_Q');
run_experiments('Profile','full','Groups',{'B3'},'Kinds',{'Q'}, ...
                'Workers',0,'Output',out);     % serial, base MATLAB
verify_run(out);
```

Run Q and G together to draw the combined Q/G B3 figure. A subset run cannot provide curves for experiments that were not simulated. The plotter reports missing/incomplete/invalid plots in `plot_audit.json` instead of substituting published data.

## Included experiments

| Group | Variable | Published curves / settings |
|---|---|---|
| A1 + A2 | LoS receive SNR −5:1:10 dB | Spatial-frequency MSE, position MSE, detection; LoS baselines and mixed NLoS JMPMD, κ=10/15 dB |
| A4 | Iterations 1–50 displayed | Position and soft-symbol MSE; QPSK/Gaussian, 0/6 dB; full cap 200 |
| B1 | Pilots 5, 10, 15, 20 | Fixed frame length 40; 6 dB; JMPMD / pilot-only / no position feedback |
| A3 | BS antennas 16:4:40 | 6 dB; JMPMD and references |
| B3 | RIS side length 8:2:16 | Fixed 25 dBm transmit power, −125 dBm noise; prior-local pilot angle initialization |

QPSK A1/A2 at 7–10 dB use 15,000 trials each; all other settings use 5,000. `experiment_catalog.m` is the executable source of truth. A5/A6, B2, B4, legacy Figure 3, and solver ablations are outside this release's main-figure scope.

## Package contents

```
run_experiments.m      resumable serial/parallel entry point
experiment_catalog.m  exact settings and deduplicated views
simulate_trial.m      one independent trial, including reference calculation
summarize_run.m        checkpoint → CSV (no failed-sample removal)
plot_results.m        published/new CSV → seven main figure groups
verify_run.m          checkpoint and completion audit
validate_release.m    small reproducibility/integration tests
src/                  six estimators, simulation, bounds, numerical helpers
private/              checkpoint, aggregation and rendering utilities
results/published/    summary/convergence CSV and 14 PDF/PNG files
metadata/             plot layout, source provenance and checksums
tests/                strict-moment tests and 60 small golden CSV samples
docs/                 protocol, math, data dictionary, figure captions
```

**No raw Monte Carlo MAT files or historical caches are distributed.** Full runs create local MAT checkpoints to allow safe resume; do not commit these generated files. The `.gitignore` excludes standard run/replot/validation locations. A full run needs substantially more disk and RAM than this small download; use an output directory with several GB available. Original runs used single-thread workers, commonly 12; elapsed time and memory depend on the host and reference recomputation.

## Reproducibility and limitations

- Trial `i` uses seed `20260927 + i`; separate random streams preserve common observations across compatible methods. Configuration IDs do not affect the seed.
- Only strict moments and the repaired local angle update are exposed. The Gaussian B3 initializer is retained; A3 uses the position prior center. There is no truth-angle initialization.
- All requested samples contribute to the ordinary mean. Numerical failures, cap hits, Taylor call diagnostics and extreme position errors remain in the output. Smoke plots are not performance evidence.
- QPSK detection uses joint-MAP SER. The single ideal-CSI curve uses the complete projected-noise covariance. Gaussian detection and A4 soft-symbol curves use MSE.
- Curves labelled BCRB retain the project's specified information conditions and **average of per-realization inverses**; they are not relabelled as a universal bound for all models.
- Exact equality of elapsed-time/memory columns is not expected. Floating-point differences across MATLAB/BLAS/platforms can change rare local solutions. A full R2023b rerun is needed to independently confirm published Monte Carlo means.

## License and attribution

Released under [MIT](LICENSE), as authorized by the code owner. Existing source attributions and the original JMPMD lineage are retained in [NOTICE.md](NOTICE.md) and `metadata/source_provenance.json`. See [CITATION.md](CITATION.md) for citation guidance. No original paper PDF or unrelated third-party algorithm distributions are bundled.
