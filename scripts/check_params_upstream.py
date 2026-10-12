"""Verify params/aerosonde.py against a local clone of byu-magicc/mavsim_public.

Usage: python scripts/check_params_upstream.py <path to mavsim_public clone>
Compares every numeric assignment that exists in both files (after the known renames).
"""
import ast
import pathlib
import sys

RENAMES = {"e": "e_oswald", "M": "M_blend"}
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from params import aerosonde as P  # noqa: E402

clone = pathlib.Path(sys.argv[1])
src = (clone / "mavsim_python/parameters/aerosonde_parameters.py").read_text()
env = {"np": __import__("numpy")}
checked = bad = 0
for node in ast.parse(src).body:
    if not (isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)):
        continue
    name = node.targets[0].id
    if name in ("e0", "e1", "e2", "e3", "e") and name != "e":
        continue
    try:
        val = eval(compile(ast.Expression(node.value), "<up>", "eval"), env)
    except Exception:
        continue
    env[name] = val
    mine = getattr(P, RENAMES.get(name, name), None)
    if mine is None or not isinstance(val, (int, float)):
        continue
    checked += 1
    if abs(mine - val) > 1e-12 * max(1, abs(val)):
        bad += 1
        print(f"MISMATCH {name}: mine={mine} upstream={val}")
print(f"checked {checked} parameters, {bad} mismatches")
sys.exit(1 if bad else 0)
