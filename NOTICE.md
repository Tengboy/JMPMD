# Authorship and source lineage

The code owner requested an MIT release retaining original author attribution. The original JMPMD authors and contributors are credited collectively in `LICENSE`; this workspace did not supply individual copyright-holder names. Existing author comments have not been replaced by invented names.

The numerical implementation derives from the original `1219_JPD_test` JMPMD experiments. The current six adapters retain their original experiment lineage in source comments. Frozen strict QPSK and Gaussian experiments were used as inputs to this release; relative source names, original hashes and release hashes are recorded in `metadata/source_provenance.json`. Those names identify historical inputs, not external runtime dependencies.

The original source tree also contained VALSE and MVALSE distributions. Their headers credit Mihai Badiu (VALSE code, March 2016; work with Thomas Hansen and Bernard Fleury) and QI Zhang / Jiang Zhu (MVALSE, based on Badiu's code). Those standalone algorithms, and the unrelated legacy estimator classes, are **not included** in this focused release. No license for those excluded distributions is asserted here.

The retained `Ainv.m` helper explicitly cites Mardia & Jupp, *Directional Statistics*, Wiley, 2000, pp. 85–86, for its inverse Bessel-ratio approximation. Its mathematical reference is retained. The repaired local angular update uses a von Mises approximation and preserves the original JMPMD update interface. Refer to `docs/METHODS.md` for implementation changes versus the original experiments.

The original paper PDF and unrelated project files are not redistributed. Paper authorship, a DOI, and a repository URL are not inferred from anonymous manuscript metadata; see `CITATION.md`.
