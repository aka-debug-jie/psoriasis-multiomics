# Third-party notices

Original project code is MIT licensed and original course prose is CC BY 4.0.
Data and bundled dependencies keep their own terms. The scope of each data grant
is recorded in [DATA_LICENSES.tsv](DATA_LICENSES.tsv); a source being downloadable
does not place it under the project's MIT licence.

The portable packages include R 4.4.2, Python, JupyterLab, R packages, and their
runtime libraries. Each package's `resources/licenses/` contains component
metadata; `component-notices.zip` preserves the full notices and recipes while
keeping the outer course paths short. Matching source archives, recipes, and patches
are offered in the same Release. Version-specific source hashes are recorded
there. The binary/runtime distribution includes GPL, LGPL, AGPL, BSD, MIT,
Apache, font licences, and Microsoft runtime redistribution terms; it is not
uniformly MIT licensed.

- R: R Core Team and contributors, GPL-2.0-or-later; retain component notices.
- SMR 1.4.3: fixed upstream commit `eedb1a62ba4c66ab0b09acd0205876955a803f18`.
  Its repository supplies a GNU GPL version 2 LICENSE file, while its startup
  banner inconsistently states “MIT License”. We retain the supplied LICENSE
  and copyright notices and provide the corresponding source and local
  modifications. SMR and modifications to it are not covered by this project's
  MIT licence. The source records the patches actually applied to each platform:
  deterministic input-index ordering for exact HEIDI |Z| ties, and Windows I/O
  and memory-reporting changes where applicable. Dependency sources are retained.
- PLINK 1.9.0: Shaun Purcell, Christopher Chang and contributors, GPL-3.0;
  fixed source commit `bceae81647e473639d6a4a99be8d23821ae3566f` corresponds to
  the September 27, 2026 release. Binary identities are separately checked.
- TwoSampleMR: Gibran Hemani and contributors, MIT. The vendored source and
  original notices are retained under `course/R/vendor/`.
- Rtsne includes software developed by Laurens van der Maaten. Its complete
  licence and attribution requirements are retained in the component notices.
- dqrng contains AGPL-licensed code; its full component licence and exact-version
  source are supplied separately from the original project code.
- Jupyter screenshots show a local JupyterLab session with original course text
  and annotation overlays. Jupyter UI copyright belongs to Project Jupyter
  contributors; JupyterLab uses the BSD-3-Clause licence. No browser address bar,
  authentication token, student work, or personal profile is included.

The original paper PDF and unreviewed supplementary attachments are not part of
this repository. Cite the paper, data providers and analytical methods used in
any derivative research, in addition to this course.

This product includes software developed by the Delft University of Technology.
