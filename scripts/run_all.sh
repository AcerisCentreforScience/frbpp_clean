#!/usr/bin/env bash
# Clean-room recompute chain for the manuscript's numbers and figure inputs.
# Run from the repository root. Set PY to the project interpreter if needed.
set -u
PY="${PY:-python}"
export PYTHONIOENCODING=utf-8
mkdir -p OUT/chain
rc=0
step() {
  name="$1"; shift
  echo "=== [$name] $* ==="
  t0=$(date +%s)
  "$@" >"OUT/chain/$name.log" 2>&1
  e=$?
  t1=$(date +%s)
  echo "--- [$name] exit=$e  $((t1 - t0))s"
  tail -3 "OUT/chain/$name.log"
  [ "$e" -ne 0 ] && rc=1
  return 0
}

START="${START:-1}"
n=0
run() { n=$((n + 1)); [ "$n" -lt "$START" ] || step "$@"; }

run 01_build_interim       "$PY" scripts/build_interim.py
run 02_event_table_audit   "$PY" scripts/event_table_dedup_audit.py
run 03_rerun_v13_dedup     "$PY" scripts/rerun_v13_event_dedup.py
run 04_gate_g4b            "$PY" scripts/gate_g4b.py
run 05_formalize_c2        "$PY" scripts/formalize_c2.py
run 06_external_v2         "$PY" scripts/c2_external_v2.py
run 07_rank_candidates     "$PY" scripts/rank_c2_candidates.py
run 08_dm_sigma_floor      "$PY" scripts/dm_sigma_floor.py
run 09_dm_floor_robust     "$PY" scripts/dm_floor_robust.py
run 10_null_battery        "$PY" scripts/c2_null_battery.py
run 11_same_beam           "$PY" scripts/same_beam_test.py --draws 400 --seed 20261011
run 12_isolation_bg        "$PY" scripts/build_isolation_data.py
run 13_known_repeater      "$PY" scripts/known_repeater_control.py
run 14_isolation_conv      "$PY" scripts/isolation_conv_565d.py
run 15_scatter             "$PY" scripts/scatter_565d.py
run 16_period              "$PY" scripts/period_565d.py
run 17_dm_z_inversion      "$PY" scripts/dm_z_inversion.py
run 18_repro_package       "$PY" scripts/repro_package.py
echo "CHAIN_RC=$rc"
exit $rc
