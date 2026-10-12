"""MATLAB/Simulink cross-validation (opt-in: needs MATLAB + Simulink; takes a few minutes).

Run with:  RUN_MATLAB_TESTS=1 python -m pytest tests/test_matlab_validation.py
It (re)generates the Python reference data if missing, runs matlab/validation/run_validation.m in batch mode
and asserts the error metrics it writes to results/phase7_validation.json.
"""

import json
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
REF = ROOT / "matlab" / "validation" / "data" / "python_reference.mat"

pytestmark = pytest.mark.skipif(
    not (os.environ.get("RUN_MATLAB_TESTS") == "1" and shutil.which("matlab")),
    reason="set RUN_MATLAB_TESTS=1 and have `matlab` on PATH to run the MATLAB/Simulink validation")


@pytest.fixture(scope="module")
def metrics():
    if not REF.exists():
        subprocess.run([sys.executable, str(ROOT / "scripts" / "make_matlab_reference.py")], check=True, cwd=ROOT)
    subprocess.run(["matlab", "-batch", "addpath('matlab/validation'); run_validation"], check=True, cwd=ROOT, timeout=1800)
    return json.loads((ROOT / "results" / "phase7_validation.json").read_text())


def test_plain_matlab_model_matches_python(metrics):
    f = metrics["functions"]
    assert f["deriv_max_abs_err"] < 1e-9
    assert f["trim_delta_max_err"] < 1e-8 and f["A_max_abs_err"] < 1e-6
    assert f["gains_max_rel_err"] < 1e-6 and f["eig_max_abs_err"] < 1e-6


@pytest.mark.parametrize("run", ["open_loop", "closed_loop_nominal", "closed_loop_wind"])
def test_simulink_time_histories_match_python(metrics, run):
    assert metrics[run]["max_abs_state_err"] < 1e-9
    assert metrics[run]["quat_norm_drift"] < 1e-5


@pytest.mark.parametrize("run", ["closed_loop_nominal", "closed_loop_wind"])
def test_autopilot_commands_and_zones_identical(metrics, run):
    assert metrics[run]["delta_cmd_max_abs_err"] < 1e-9
    assert metrics[run]["zone_mismatches"] == 0
