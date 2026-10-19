# Reproduction protocol

## Supported environment

- Tested: MATLAB R2023b, Windows, one computational thread per worker.
- Required for serial computation: base MATLAB with its normal JVM. No Optimization, Statistics, Communications or external MATLAB packages are needed.
- Optional: Parallel Computing Toolbox for `Workers > 0`. The runner uses a local process pool. Default `Workers=0` executes sequentially.
- Plotting: a Chinese-capable font (SimSun, Noto Serif CJK SC, Source Han Serif SC, Microsoft YaHei, PingFang SC or a user-specified font). Numeric tick labels use Times New Roman or MATLAB's substitution.
- Python is not needed to simulate or plot. `tests/check_distribution.py` optionally verifies release files using Python's standard library.
- Linux/macOS are intended to be portable but were not tested for numerical equivalence in this release. The Windows-only `memory` diagnostic becomes NaN if unavailable; this does not fail a trial.

Use the release's directory on the MATLAB path, not `genpath` of a workspace containing multiple older versions. `setup_jmpmd` adds only this release and its numerical helpers. Restart MATLAB or `restoredefaultpath` first if other versions of identically named helper functions are already loaded.

## Catalog and counts

`[conditions,views]=experiment_catalog` returns the full catalog without launching anything. `conditions` identifies independent estimator/configuration runs. `views` maps the same results to experiment figures. Deduplication is performed on the complete configuration, including initialization.

| Selection | Unique configs | Trials | Views |
|---|---:|---:|---:|
| QPSK, all main groups | 110 | 750000 | 116 |
| Gaussian, all main groups | 110 | 550000 | 116 |
| Combined | 220 | 1300000 | 232 |

A4 shares A1's LoS JMPMD 0/6 dB observations. A3 at BS=32 and B1 at Tp=10 share A1's 6 dB results when run together. Separate invocations into separate output folders necessarily recompute overlapping settings; run all groups together for global deduplication. A group selection only computes references requested by that group.

`smoke` replaces each target by 2 trials and the 200-iteration cap by 4. It exercises all group/method paths with 440 estimator trials, 88 reference evaluations and 8 full short trajectories. It is a pipeline test only.

## Randomization

All trials use `seed = 20260927 + trial_id`. The condition number never changes observations. The scene uses `mt19937ar` streams:

| Offset | Draws |
|---:|---|
| 0 | User position |
| 2000000 | Symbols (QPSK or circular complex Gaussian) |
| 4000000 | RIS probes |
| 6000000 | Additive noise |
| 8000000 | NLoS master template |
| 10000000 | Estimator RNG initialization |

Methods at the same physical setting share observations; Q and G have distinct symbol laws. Different array sizes change draw shapes, so pairing should not be interpreted as identical observations across different geometries. NLoS branches use a common master template so branch selection does not shift other random draws.

## Checkpoints and resume

```matlab
out=fullfile(pwd,'runs','full');
run_experiments('Profile','full','Workers',4,'MaxBatches',10,'Output',out);
% Run the identical selection again; worker count/batch size may be changed.
run_experiments('Profile','full','Workers',4,'Output',out);
```

The current output path is locked before initializing or modifying run metadata. `process.json` records the actual MATLAB PID and host. Existing committed `chunks/cNNN_first_last.mat` are never overwritten. Source signature, configuration, MATLAB version, trial order, seed, record dimensions and duplicate/gap conditions are checked before continuation. A crash can require redoing an uncommitted batch, not already committed samples.

For graceful stop, create an empty `STOP_AFTER_CHECKPOINT` file inside the output folder. Remove it before continuing. A `MaxBatches` stop automatically removes its lock. An abrupt power loss may leave `RUNNING.lock`: verify on the recorded host that the recorded MATLAB process is no longer alive before deleting **that lock only**. Never remove a live run's lock or kill unrelated MATLAB processes. Failed-sample records do not justify deleting checkpoints or silently rerunning selected samples.

`manifest.mat` contains configs/metadata; `chunks/*.mat` contain per-trial results and selected traces. `status.json` is progress metadata and can be regenerated from chunks. Raw MAT files are local generated outputs, excluded from the release. The source signature covers executable MATLAB files but not documentation or the published CSV.

Use a fresh output directory for changed source/config/profile/MATLAB. The directory name is not inferred from group selections: specify `Output` whenever starting an independent subset run. The runner refuses to write inside `results/published`.

## Export and audit

`summarize_run(out)` can rebuild CSV without simulation. In incomplete runs `Pending>0`, and primary means for unfinished views are NaN. `verify_run(out)` requires all trials committed and checks structure; its JSON separately reports numerical failures. A structurally complete run with failed samples is **not** a claim of numerical success.

`plot_results(out,figDir)` reads only that directory's CSV. Incomplete groups, duplicate points and nonfinite/nonpositive logarithmic values are reported in `plot_audit.json`; no replacement by bundled results, tail removal, pseudo-count or zero-to-positive conversion is performed. Figures whose required curves are all available are exported. `results/published` retains the verified reference figures regardless of what a new short run can plot.

Elapsed seconds and memory in the historical CSV are measurements from the original machines/runs; the release recomputes references and may use different batching. They will not match new timing values. Published mean/SEM values came from complete formal experiments; the provided release checks are smaller regression/integration checks, not a repeat of all 1.3 million trials.
