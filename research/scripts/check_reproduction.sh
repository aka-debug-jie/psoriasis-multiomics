#!/usr/bin/env bash
# Repeatable acceptance checks for the analyses that have actually run.
set -euo pipefail
cd "$(dirname "$0")/../.."
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
bash research/scripts/bootstrap_env.sh
bash research/scripts/bootstrap_data_env.sh
bash research/scripts/bootstrap_sc_env.sh
research/envs/data/bin/python -m unittest discover -s research/tests -p 'test_*.py'
research/envs/r442/bin/Rscript --vanilla research/tests/test_mr_core.R
research/envs/r442/bin/Rscript --vanilla research/tests/test_mr_available.R
research/envs/r442/bin/Rscript --vanilla research/tests/test_author_scores.R
research/envs/r442_sc/bin/Rscript --vanilla research/tests/test_sc_de_math.R
research/envs/data/bin/python research/tests/validate_results.py
research/envs/data/bin/python research/tests/validate_completion.py
research/envs/data/bin/python research/tests/validate_sc_embedding.py
research/envs/data/bin/python research/tests/validate_sc_de.py
research/envs/data/bin/python research/tests/validate_new_analyses.py
if [[ -f research/release/figures_manifest.json ]]; then
  research/envs/data/bin/python research/tests/validate_release.py
fi
if [[ -f research/results/covariance_resolution/run.json ]]; then
  research/envs/data/bin/python research/tests/validate_covariance_resolution.py
fi
if [[ -f research/final/run.json ]]; then
  research/envs/data/bin/python research/tests/validate_final_research.py
fi
