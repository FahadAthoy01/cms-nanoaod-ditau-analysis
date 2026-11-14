import numpy as np
import ROOT
import src.reconstruct_ditau  # Registers C++ helper functions in ROOT interpreter


def test_transverse_mass_back_to_back():
    mt = ROOT.compute_mt(40.0, 0.0, 40.0, float(np.pi))
    assert np.isclose(mt, 80.0, atol=1e-3)


def test_transverse_mass_collinear():
    mt = ROOT.compute_mt(35.0, 1.2, 50.0, 1.2)
    assert np.isclose(mt, 0.0, atol=1e-4)


def test_delta_r_orthogonal():
    dr = ROOT.compute_delta_r(0.0, 0.0, 0.3, 0.4)
    assert np.isclose(dr, 0.5, atol=1e-4)


def test_visible_mass_on_shell():
    m_vis = ROOT.compute_vis_mass(45.6, 0.0, 0.0, 0.10566, 45.6, 0.0, float(np.pi), 1.77686)
    assert np.isclose(m_vis, 91.2, atol=0.2)
