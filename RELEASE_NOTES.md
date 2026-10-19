# 2026.09.30 — main-figure reproduction release

- Includes the seven selected thesis figure groups: A1/A2 theta/position/detection, A4 convergence, B1 pilots, A3 BS antennas and B3 RIS side length, with strict QPSK and Gaussian implementations.
- Removes dependence on historical comparison/reference MAT caches. Fresh experiments calculate their own observations, bounds and ideal-CSI results.
- Adds a deduplicated experiment catalog, serial/optional parallel runner, immutable batches, restart checks, failure-preserving CSV export and one plotting entry point.
- Bundles only lightweight published CSV and PDF/PNG, plus small regression fixtures. No full Monte Carlo result MAT files are distributed.
- Preserves numerical update order; Gaussian estimators include the optional B3 angle-initialization hook. Removes only an unused Q-P legacy object constructor and an unselected legacy angle-solver dispatch.
- Adds portable Chinese font detection, bilingual quick starts, model/metric and implementation notes, provenance hashes, MIT license and retained lineage attribution.

See `docs/VALIDATION.md` for measured verification and its limits. The release prepares local files for publication; repository hosting, a DOI and complete personal/bibliographic attribution are maintainer metadata, not inferred values.
