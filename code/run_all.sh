#!/bin/sh
# V0.2: reproduce the complete study from this package root. Each step logs to logs/<step>.log and
# records its exit code; the analysis still runs when the main run reports unresolved records, so
# incomplete summaries and the failure ledger remain available. Final exit status is nonzero if any
# step failed or any record is unresolved.
# Requires: Python 3, cypari2 (PARI/GP; recorded run used 2.19.1), mpmath, numpy, matplotlib.
mkdir -p logs tests results
rc_all=0
step() { name=$1; shift; "$@" > logs/$name.log 2>&1; rc=$?; echo "$name exit=$rc" | tee -a logs/exit_codes.txt; [ $rc -ne 0 ] && rc_all=1; return 0; }
: > logs/exit_codes.txt
step selftests          python3 code/selftests.py inputs/allgens inputs/handoff_reference/reduced_basis_shape_extract.csv tests/selftests.json
step negative_controls  python3 code/negative_controls.py tests/negative_controls.json
step integration_tests  python3 code/integration_tests.py inputs/allgens tests/integration_work tests/integration_tests.json
step run_study          python3 code/run_study.py inputs/allgens inputs/handoff_reference/source_counts.json results --workers 2
step retry_exercise     python3 code/retry_path_exercise.py results/per_curve_results.jsonl tests/retry_path_exercise.json
step analyze            python3 code/analyze.py results
step plots              python3 code/plots.py results
echo "run_all exit=$rc_all" | tee -a logs/exit_codes.txt
exit $rc_all
