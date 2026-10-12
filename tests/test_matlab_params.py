"""The MATLAB parameter file is generated from params/aerosonde.py; keep them identical (runs without MATLAB)."""

import pathlib
import re

import pytest

from params import aerosonde as P

M_FILE = pathlib.Path(__file__).resolve().parents[1] / "matlab" / "uav_params.m"


def _parse():
    vals = {}
    for line in M_FILE.read_text(encoding="utf-8").splitlines():
        m = re.fullmatch(r"P\.(\w+) = ([^;]+);", line.strip())
        if m:
            vals[m.group(1)] = float(m.group(2))
    return vals


def test_matlab_params_match_python_exactly():
    vals = _parse()
    py = {n: float(getattr(P, n)) for n in dir(P) if not n.startswith("_") and isinstance(getattr(P, n), (int, float))}
    assert set(vals) == set(py), set(vals) ^ set(py)
    for k, v in py.items():
        assert vals[k] == v, k          # repr() round trip: bit-identical


def test_matlab_params_file_is_marked_generated():
    assert "GENERATED" in M_FILE.read_text(encoding="utf-8")
