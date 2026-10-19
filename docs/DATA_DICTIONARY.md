# Lightweight result tables

## summary.csv

One row is one experiment view of a condition, not an additional independent run. The published table has 232 rows. Match curves by `(Group, Kind, Series, X)`. `Condition` is local to its source dataset; published rows from different historical runs must not be joined by that ID alone.

| Columns | Meaning |
|---|---|
| Group / Kind / Method / Series | A1/A2/A3/A4/B1/B3; Q/G; J/P/W; series J/P/W or K10/K15 |
| X | Sweep value: SNR, BS antennas, pilots, or RIS side length |
| SNR / NoiseMode / ActualSNR | Config SNR; `snr` or `fixed`; measured receive LoS SNR. For fixed mode SNR is an unused config value. |
| NB / M / NR / Tp / Td | BS antennas, RIS count, RIS elements, pilot/data slots |
| NLoS / Kappa / InputKappa / ProjectedKappa | Model and preprojection target/measured ratio; projected ratio |
| Strict / InitMode / AngleSolver / GainPilotInit | Estimator settings; all published rows are strict/local |
| Target / N / Pending | Requested, committed and uncommitted counts |
| Theta / Position / Data / ThetaPost | Main sample means; units and modulation meaning in METHODS.md |
| ThetaSEM / PositionSEM / DataSEM | Ordinary sample standard error of the mean |
| PositionBound / ThetaBound / DataBound | Averaged reference values; structurally unavailable values are NaN |
| CSIFullSER | Published QPSK ideal-CSI full-covariance joint MAP SER |
| CSISER | Historical diagonal-noise diagnostic; not plotted as a second benchmark |
| Failed / PrimaryNonfinite | Whole trial failures; nonfinite primary scalar entries (up to 3 per trial) |
| AtCap / Unconverged | Iteration cap hits; trials not satisfying convergence. A converged trial may also hit the cap. |
| TaylorCalls / TaylorFailedCalls / TaylorAffectedTrials | Call counts and affected-trial counts; not synonymous with whole-trial failure |
| LargePositionTrials / LargePositionUsers | At least one user, or user instances, with squared position error >1 m² |
| UnknownPositionUsers | Nonfinite final user position errors |
| ReferenceFailed / ReferenceMissing | Numerical reference failures / missing required references |
| Iterations / EstimatorSeconds / TotalSeconds / BoundSeconds | Mean iteration count and elapsed time diagnostics |
| MaxBoundResidual / MaxWorkerMemoryGB | Reference inverse consistency diagnostic and sampled worker memory |
| SymbolErrors / SymbolDecisions | QPSK aggregate count; not applicable to Gaussian |
| SourceDataset / SourceRun / SourceCondition / SourceSignature | Historical provenance where recorded; portable identifiers only |

References are requested for LoS J conditions in A1/A3; shared conditions may expose the same values in B1/A4 views. NaN in non-applicable reference columns is not a primary estimator failure. When aggregating failures or total work, deduplicate by the physical condition; do not simply sum across shared views.

## convergence.csv

800 published rows: Q/G × 0/6 dB × 200 iterations, 5000 trajectories per setting. Columns: `Kind, SNR, Iteration, N, G_MSE, G_NMSE, X_MSE, Theta, Position, G_NMSE_SEM, X_MSE_SEM, Position_SEM, Active` plus provenance. `X_MSE` is the soft posterior-mean error for both modulations. `Active` counts trajectories performing that iteration. Means and SEMs use all trajectories with the final value held after stopping, not only the active subset.

## Local generated checkpoints

`manifest.mat` and `chunks/*.mat` are created only by a new run. Record columns are listed in `private/record_columns.m` and the JSON manifest. Each chunk includes trial IDs, seeds, source signature, 43 scalar diagnostics, optional full traces, warning identifiers and exception text. Timing/memory entries in `tests/golden_samples.csv` are historical diagnostics and are explicitly excluded from numerical regression comparisons.

## Published snapshots versus fresh runs

The CSV/PDF/PNG under `results/published` are frozen lightweight outputs from complete formal experiments. They are not inputs to `simulate_trial`. The plot layout JSON contains selectors and styling, not substitute curve y-values. `plot_results` resolves every point from the chosen CSV directory. Missing or invalid results are reported, not silently replaced.
