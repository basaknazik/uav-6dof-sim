"""Visualisation: 3D trajectory plot and an attitude animation (matplotlib only).

Plot axes are East (x), North (y), Up (z) so the picture looks like a map; all state data are NED
(position [pn, pe, pd]); altitude = -pd. The aircraft is a simple line drawing in the body frame
(x nose, y right wing, z down) rotated into NED with R_b_to_v = R_v_to_b^T.
"""

import numpy as np

from sim import rotations as rot

# Simple airframe in body axes, unit = fraction of `scale` metres: fuselage, wing, horizontal tail, fin.
_BODY_LINES = [
    np.array([[1.0, 0, 0], [-1.0, 0, 0]]),                       # fuselage (nose at +x)
    np.array([[0.1, -1.6, 0], [0.1, 1.6, 0]]),                   # main wing
    np.array([[-0.9, -0.6, 0], [-0.9, 0.6, 0]]),                 # horizontal tail
    np.array([[-0.9, 0, 0], [-0.9, 0, -0.5]]),                   # vertical fin (up = -z body)
]


def aircraft_lines(euler, pos_ned=(0.0, 0.0, 0.0), scale=1.0, z_scale=1.0):
    """Line segments of the airframe drawing in plot axes (E, N, Up).

    euler = (phi, theta, psi) [rad]; pos_ned [m] (CG); scale [m] multiplies the unit drawing
    (fuselage is 2*scale long). z_scale multiplies the glyph's vertical extent about its CG: use 1/zf when the
    plot's vertical axis is stretched by a factor zf, so the airframe looks undistorted. Returns a list of (2, 3) arrays.
    """
    R_b_to_v = rot.euler_to_R(*euler).T
    out = []
    for seg in _BODY_LINES:
        rel = (R_b_to_v @ (scale * seg).T).T                 # offsets from the CG in NED
        pos = np.asarray(pos_ned, float)
        out.append(np.column_stack([pos[1] + rel[:, 1], pos[0] + rel[:, 0], -pos[2] - z_scale * rel[:, 2]]))   # (E, N, Up)
    return out


def euler_from_quat_series(X):
    """(N, 3) array of (phi, theta, psi) from the quaternion columns 6:10 of a state history."""
    return np.array([rot.quat_to_euler(q) for q in X[:, 6:10]])


def plot_trajectory_3d(ax, X, n_aircraft=12, scale=30.0, z_exaggeration=8.0, label="trajectory"):
    """Draw a 3D trajectory (E, N, altitude) with airframe glyphs along it. X: state history (N, >=10).

    The vertical axis is stretched by `z_exaggeration` (altitude changes are tiny next to the ground track);
    glyphs are compensated so the airframe is not distorted.
    """
    E, N, H = X[:, 1], X[:, 0], -X[:, 2]
    ax.plot(E, N, H, lw=1.2, label=label)
    ax.plot(E, N, np.full_like(H, H.min() - 10), lw=0.6, color="gray", alpha=0.6)   # ground shadow
    eul = euler_from_quat_series(X)
    for k in np.linspace(0, len(X) - 1, n_aircraft).astype(int):
        for seg in aircraft_lines(eul[k], X[k, 0:3], scale, 1.0 / z_exaggeration):
            ax.plot(*seg.T, color="k", lw=1.0)
    ax.scatter(E[0], N[0], H[0], color="g", s=25)
    ax.scatter(E[-1], N[-1], H[-1], color="r", s=25)
    ax.set(xlabel="East [m]", ylabel="North [m]", zlabel="Altitude [m]")
    span = max(np.ptp(E), np.ptp(N), 1.0)
    ax.set_box_aspect((max(np.ptp(E), 0.05 * span), max(np.ptp(N), 0.05 * span), z_exaggeration * (np.ptp(H) + 20)))
    return ax


def make_animation(X, t, path, fps=20, stride=10, scale=40.0, dpi=70):
    """Two-panel GIF: 3D flight path with the airframe, and a large attitude view (nose/wings) with Euler angles.

    X: state history (N, >=13 quaternion model); t: time [s]; path: output .gif; stride: frames skipped.
    Returns the number of frames written.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    eul = euler_from_quat_series(X)
    frames = np.arange(0, len(X), stride)
    fig = plt.figure(figsize=(9, 4.2))
    a1 = fig.add_subplot(1, 2, 1, projection="3d")
    a2 = fig.add_subplot(1, 2, 2, projection="3d")
    E, N, H = X[:, 1], X[:, 0], -X[:, 2]
    a1.plot(E, N, H, color="lightgray", lw=0.8)
    a1.set(xlabel="E [m]", ylabel="N [m]", zlabel="h [m]")
    a1.set_xlim(E.min() - 100, E.max() + 100)
    a1.set_ylim(N.min() - 100, N.max() + 100)
    a1.set_zlim(H.min() - 20, H.max() + 20)
    zf = 6.0
    a1.set_box_aspect((np.ptp(E) + 200, np.ptp(N) + 200, zf * (np.ptp(H) + 40)))
    trail, = a1.plot([], [], [], color="C0", lw=1.5)
    glyph1 = [a1.plot([], [], [], color="k", lw=1.5)[0] for _ in _BODY_LINES]
    glyph2 = [a2.plot([], [], [], color="k", lw=2.5)[0] for _ in _BODY_LINES]
    for ax in (a2,):
        ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_zlim(-2, 2)
        ax.set(xlabel="East", ylabel="North", zlabel="Up", title="Attitude (nose = long axis)")
        ax.set_box_aspect((1, 1, 1))
    txt = fig.text(0.5, 0.02, "", ha="center")

    def update(k):
        i = frames[k]
        trail.set_data_3d(E[:i + 1], N[:i + 1], H[:i + 1])
        for g, seg in zip(glyph1, aircraft_lines(eul[i], X[i, 0:3], scale, 1.0 / zf)):
            g.set_data_3d(*seg.T)
        for g, seg in zip(glyph2, aircraft_lines(eul[i], (0, 0, 0), 1.0)):
            g.set_data_3d(*seg.T)
        phi, theta, psi = np.rad2deg(eul[i])
        txt.set_text(f"t = {t[i]:6.1f} s   roll {phi:6.1f}°  pitch {theta:5.1f}°  yaw {psi:6.1f}°   alt {H[i]:6.1f} m")
        return [trail, *glyph1, *glyph2, txt]

    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    anim.save(path, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return len(frames)
