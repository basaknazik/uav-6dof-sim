"""Phase 6: the visualisation helpers draw the airframe with the right geometry (this is where sign/frame bugs hide)."""

import numpy as np
import pytest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sim import rotations as rot
from sim import visualize as vz

D = np.deg2rad
NOSE, RIGHT_TIP, FIN_TOP = 0, 1, 3      # indices into the segment list; point 0 of fuselage = nose (+x body)


def _nose(lines):
    return lines[0][0]


def _right_wingtip(lines):
    return lines[1][1]


def test_level_northbound_airframe_geometry():
    L = vz.aircraft_lines((0, 0, 0), scale=10.0)                     # plot axes (E, N, Up)
    np.testing.assert_allclose(_nose(L), [0, 10, 0], atol=1e-12)     # nose points North
    np.testing.assert_allclose(_right_wingtip(L), [16, 1.0, 0], atol=1e-12)   # right wing points East (+E)
    np.testing.assert_allclose(L[3][1][2], 5.0, atol=1e-12)          # fin tip is ABOVE the CG (Up)


def test_yaw_90_points_nose_east_and_right_wing_south():
    L = vz.aircraft_lines((0, 0, D(90)), scale=10.0)
    np.testing.assert_allclose(_nose(L), [10, 0, 0], atol=1e-12)
    assert _right_wingtip(L)[1] < 0                                  # right wing now points South


def test_pitch_up_raises_nose_and_roll_right_drops_right_wing():
    assert _nose(vz.aircraft_lines((0, D(20), 0), scale=10.0))[2] == pytest.approx(10 * np.sin(D(20)))
    tip = _right_wingtip(vz.aircraft_lines((D(90), 0, 0), scale=10.0))
    assert tip[2] == pytest.approx(-16.0)                            # right wing down
    assert tip[0] == pytest.approx(0.0, abs=1e-12)


def test_position_and_altitude_convention():
    L = vz.aircraft_lines((0, 0, 0), pos_ned=(100.0, 200.0, -50.0), scale=1.0)   # N=100, E=200, h=50
    cg = 0.5 * (L[0][0] + L[0][1])
    np.testing.assert_allclose(cg, [200.0, 100.0, 50.0], atol=1e-12)


def test_z_scale_only_stretches_vertical_extent_about_the_cg():
    a = vz.aircraft_lines((D(10), D(20), D(30)), (0, 0, -100.0), 5.0, 1.0)
    b = vz.aircraft_lines((D(10), D(20), D(30)), (0, 0, -100.0), 5.0, 0.25)
    for sa, sb in zip(a, b):
        np.testing.assert_allclose(sb[:, :2], sa[:, :2])
        np.testing.assert_allclose(sb[:, 2] - 100.0, 0.25 * (sa[:, 2] - 100.0), atol=1e-12)


def test_lines_match_rotation_matrix_for_random_attitude():
    rng = np.random.default_rng(0)
    for _ in range(20):
        e = (rng.uniform(-1, 1), rng.uniform(-1.2, 1.2), rng.uniform(-3, 3))
        L = vz.aircraft_lines(e, scale=1.0)
        v_ned = rot.euler_to_R(*e).T @ np.array([1.0, 0, 0])         # body x-axis expressed in NED
        np.testing.assert_allclose(_nose(L), [v_ned[1], v_ned[0], -v_ned[2]], atol=1e-12)


def test_euler_from_quat_series_roundtrip():
    eul = np.array([[0.1, 0.2, 0.3], [-0.4, 0.5, 2.0]])
    X = np.zeros((2, 13))
    X[:, 6:10] = [rot.euler_to_quat(*e) for e in eul]
    np.testing.assert_allclose(vz.euler_from_quat_series(X), eul, atol=1e-12)


def test_trajectory_plot_and_animation_write_files(tmp_path):
    n = 60
    X = np.zeros((n, 13))
    X[:, 1] = np.linspace(0, 600, n)                                   # flying East
    X[:, 2] = -100.0 - np.linspace(0, 20, n)                           # climbing
    X[:, 6:10] = [rot.euler_to_quat(0.3 * np.sin(k / 8), 0.05, D(90)) for k in range(n)]
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    vz.plot_trajectory_3d(ax, X, n_aircraft=5)
    assert len(ax.lines) >= 1 + 1 + 5 * 4                              # path + shadow + 5 airframes x 4 segments
    fig.savefig(tmp_path / "t.png"); plt.close(fig)
    assert (tmp_path / "t.png").stat().st_size > 1000
    frames = vz.make_animation(X, np.arange(n) * 0.5, str(tmp_path / "a.gif"), fps=5, stride=10, dpi=40)
    assert frames == 6
    assert (tmp_path / "a.gif").stat().st_size > 1000
    assert (tmp_path / "a.gif").read_bytes()[:6] in (b"GIF87a", b"GIF89a")
