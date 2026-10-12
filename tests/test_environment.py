"""Phase 0 smoke tests: the toolchain is importable and the package layout works.

These are deliberately not physics tests; physics verification starts in Phase 1.
"""

import numpy as np
import scipy
import matplotlib


def test_scientific_stack_imports():
    assert int(np.__version__.split(".")[0]) >= 2
    assert int(scipy.__version__.split(".")[0]) >= 1
    assert matplotlib.__version__


def test_project_packages_importable():
    import sim
    import params.aerosonde  # noqa: F401

    assert "NED" in sim.__doc__


def test_matplotlib_headless_backend():
    """Plots must render without a display (CI / batch runs)."""
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    plt.close(fig)
