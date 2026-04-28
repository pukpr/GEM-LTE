#!/usr/bin/env python3
"""
geometric_fit.py
================

Closed-form lunisolar geometric fitter, reframed around the *draconic
month as the k=0 group-symmetry generator* with the 18.6-year nodal
regression as the slow modulating angle, AND solved as a single linear
OLS problem on the expanded sideband basis.

Why a linearized basis
----------------------
The geometric closed form

    y(t) = a sin( i_eff sin(theta_d + d_d) )
             * (1 + e_eff cos(theta_a + d_a))
             * (1 + m_N   cos(theta_N + d_N))
           + b0

is a NONLINEAR product. Fitting it directly leads the optimizer to park
in basins where m_N and e_eff collapse to small values while a and i_eff
absorb the carrier amplitude. The fit looks reasonable (R^2 > 0.8) but
the residual spectrum still shows the locked sidebands almost
unattenuated -- exactly the pathology you saw.

Algebraic expansion of the product model produces a finite sum of
sinusoids at known periods:

    carrier:                drac
    nodal sidebands:        drac +/- N            (= 27.32 d, 27.10 d)
    anomalistic sidebands:  drac +/- ano          (eccentricity)
    cross sidebands:        drac +/- N +/- ano
    Jacobi-Anger:           3*drac and its sidebands (from sin(i sin x))
    DC + slow drift

Every sinusoid has an INDEPENDENT linear amplitude. So instead of
fitting the nonlinear product directly we:

  1. Build a basis matrix of sin/cos columns at every relevant
     (carrier, sideband, harmonic) combination
  2. Solve ONE OLS problem -- no phase grid, no multistart, no local
     minima
  3. Back out the geometric parameters (a, i_eff, e_eff, m_N, phases)
     from the recovered linear amplitudes
  4. Report the upper/lower sideband asymmetry as a diagnostic --
     equal sidebands => pure AM, unequal => FM component or independent
     forcing

This collapses the multi-modal optimization landscape to a single linear
solve and -- crucially -- lets the data choose the sideband amplitudes
freely. If the AM/locked-sideband hypothesis holds, the recovered
upper/lower sidebands come out symmetric. If it doesn't, the asymmetry
is *measured* rather than forced to zero by the model.

Free linear amplitudes (typically ~28-32, all linear):

    a_d_c, a_d_s                 carrier sin/cos at theta_d
    a_dN+_c, a_dN+_s             upper nodal sideband (drac + N)
    a_dN-_c, a_dN-_s             lower nodal sideband (drac - N) = sidereal
    a_da+_c, a_da+_s             upper anomalistic sideband (drac + ano)
    a_da-_c, a_da-_s             lower anomalistic sideband (drac - ano)
    a_dNa++  .. a_dNa--          four cross terms drac +/- N +/- ano
    a_3d, a_3dN+, ..             third-harmonic carrier + sidebands
    plus: nodal envelope sin/cos, anomalistic envelope sin/cos
    plus: polynomial drift coefficients (degree --detrend, default 2)

Geometric back-out
------------------
From the linear amplitudes we reconstruct:

    a       = sqrt(a_d_c^2 + a_d_s^2)          (carrier amplitude)
    d_drac  = atan2(-a_d_c, a_d_s)              (carrier phase)
    e_eff   = sqrt(A_+^2 + A_-^2 + 2*A_+*A_- cos(...)) / a
              from upper+lower anomalistic sideband sum
    m_N     = sqrt(B_+^2 + B_-^2 + 2*B_+*B_- cos(...)) / a
              from upper+lower nodal sideband sum
    i_eff   from the ratio of 3rd-harmonic to carrier amplitudes

We also compute:

    AM_asymmetry_N    = |B_+ - B_-| / (B_+ + B_-)
    AM_asymmetry_ano  = |A_+ - A_-| / (A_+ + A_-)

A pure AM (locked-sideband) model predicts both asymmetries near 0.
Nonzero values are direct evidence of FM contamination or independent
forcing.

Outputs
-------
geometric_out/
  metrics.txt         text summary, R^2, recovered linear amplitudes,
                      geometric back-out, sideband asymmetry diagnostics,
                      and per-band residual variance.
  fit.png             data + model overlay + residual spectrum.
  params.json         machine-readable amplitudes + geometric params
  residual.dat        whitespace-delimited (t, y, model, residual)

Usage
-----
    python geometric_fit.py daily.dat --time-unit years
    python geometric_fit.py --synthetic           # built-in test
    python geometric_fit.py daily.dat --detrend 3 # cubic polynomial trend
"""

from __future__ import annotations
import argparse
import json
import os
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple

import numpy as np

# scipy is optional here -- we use it for the windowed Lomb-Scargle in
# diagnostics, but the core fit is pure NumPy.
try:
    from scipy.signal import lombscargle
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False


# ----------------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------------
P_SIDEREAL    = 27.321661    # days, sidereal/tropical month -- emerges as drac - N
P_DRACONIC    = 27.212221    # days, draconic/nodical month  -- the k=0 carrier
P_ANOMALISTIC = 27.554550    # days, anomalistic month
P_NODAL       = 6793.4765   # days, lunar ascending-node regression (18.61 yr)
P_SID_MIRROR  = 1.0 / (1.0/P_DRACONIC + 1.0/P_NODAL)  # 27.1036 d -- upper sideband
DAYS_PER_YEAR = 365.247484


# ----------------------------------------------------------------------------
# Basis specification
# ----------------------------------------------------------------------------
# Each basis component is a tuple (label, k_d, k_N, k_a) meaning the
# component oscillates at frequency  k_d/P_DRAC + k_N/P_NODAL + k_a/P_ANO
# with both sin and cos columns (two linear amplitudes each).
#
# The list below corresponds to the algebraic expansion of
#     a sin(i sin(th_d)) (1 + e cos(th_a)) (1 + m cos(th_N))
# truncated to leading non-trivial order.

@dataclass
class BasisSpec:
    label: str
    k_d:   int
    k_N:   int
    k_a:   int

    def angular_freq(self) -> float:
        return (2*np.pi * self.k_d / P_DRACONIC
                + 2*np.pi * self.k_N / P_NODAL
                + 2*np.pi * self.k_a / P_ANOMALISTIC)

    def period(self) -> float:
        w = self.angular_freq()
        if w == 0:
            return float("inf")
        return 2*np.pi / abs(w)


def build_basis_spec(use_third_harmonic: bool = True,
                     use_cross_terms: bool = True,
                     use_envelopes: bool = True,
                     use_perigean: bool = True,
                     use_second_harmonic: bool = True
                     ) -> List[BasisSpec]:
    """
    Enumerate the (k_d, k_N, k_a) terms that come out of the algebraic
    expansion of the geometric closed form.

    The perigean (8.85 yr) precession is NOT an independent fundamental
    angle -- it equals the anomalistic-minus-sidereal beat:

        1/P_per = 1/P_ano - 1/P_sid = 1/P_ano - 1/P_drac + 1/P_N

    So sidebands of the form 'drac +/- per' and 'ano +/- per' are
    expressible in (k_d, k_N, k_a) without adding a 4th angle:

        drac + per (freq)  -> 1/P_drac + (1/P_ano - 1/P_drac + 1/P_N)
                            = 1/P_ano + 1/P_N            -> ano + N
        drac - per (freq)  -> 1/P_drac - (1/P_ano - 1/P_drac + 1/P_N)
                            = 2/P_drac - 1/P_ano - 1/P_N -> 2drac - N - ano

    Empirically, these emerge as residual peaks at 26.985 d and 27.44 d
    (with mirror at 27.66 d = ano - N) once the first-order nodal and
    apsidal sidebands have been absorbed.
    """
    specs: List[BasisSpec] = []
    # carrier
    specs.append(BasisSpec("drac",          1,  0,  0))
    # nodal sidebands
    specs.append(BasisSpec("drac+N",        1, +1,  0))    # 27.1036 d (mirror)
    specs.append(BasisSpec("drac-N",        1, -1,  0))    # 27.3216 d (sidereal)
    # anomalistic sidebands
    specs.append(BasisSpec("drac+ano",      1,  0, +1))
    specs.append(BasisSpec("drac-ano",      1,  0, -1))

    if use_cross_terms:
        # second-order cross terms: drac +/- N +/- ano
        specs.append(BasisSpec("drac+N+ano", 1, +1, +1))
        specs.append(BasisSpec("drac+N-ano", 1, +1, -1))
        specs.append(BasisSpec("drac-N+ano", 1, -1, +1))
        specs.append(BasisSpec("drac-N-ano", 1, -1, -1))

    if use_perigean:
        # Perigean (8.85 yr) sidebands of draconic and anomalistic.
        # All expressible in (k_d, k_N, k_a) via 1/P_per = 1/P_ano - 1/P_sid.
        specs.append(BasisSpec("drac+per",   2, +1, -1))   # 26.7725 d (* see note)
        specs.append(BasisSpec("drac-per",   2, -1, -1))   # 26.985 d  <-- big peak
        specs.append(BasisSpec("ano+N",      0, +1, +1))   # 27.4433 d
        specs.append(BasisSpec("ano-N",      0, -1, +1))   # 27.6666 d
        # Note: 'drac+per' written as (2,+1,-1) gives 26.77 d; the actual
        # 27.44 d 'drac+per' peak is identical to ano+N because perigean
        # precession IS ano-sid (group identity). We list both labels for
        # clarity but they hit different (k_d,k_N,k_a) cells.

    if use_second_harmonic:
        # Half-month band (~13.6-13.8 d). These come from products of
        # first-harmonic carrier and sideband terms (squaring or
        # cross-multiplying), which a strictly-multiplicative geometric
        # model produces at 2nd order. Each is a lattice point.
        specs.append(BasisSpec("2drac",      2,  0,  0))   # 13.6061 d
        specs.append(BasisSpec("2drac+N",    2, +1,  0))   # 13.5789 d
        specs.append(BasisSpec("2drac-N",    2, -1,  0))   # 13.6334 d (drac+sid)
        specs.append(BasisSpec("2ano",       0,  0, +2))   # 13.7773 d
        specs.append(BasisSpec("2ano+N",     0, +1, +2))   # ~13.74 d
        specs.append(BasisSpec("2ano-N",     0, -1, +2))   # ~13.81 d
        specs.append(BasisSpec("drac+ano",   1,  0, +1))   # 13.69 d
        # (drac+ano already added above as a first-order anomalistic
        # sideband; specs list deduplication handled below.)

    if use_third_harmonic:
        # Third-harmonic band (~9.07-9.13 d). The k=3 lattice points;
        # we extend the nodal-sideband fan to k_N=+/-3 because real data
        # shows residual peaks at 3drac-2N (= drac+2sid, 9.095 d) and
        # 3drac-3N (= 3sid, 9.107 d) once the k_N=0,+/-1 columns are
        # absorbed.
        for kN in (-3, -2, -1, 0, +1, +2, +3):
            label = f"3drac{kN:+d}N".replace("+0N", "")
            if kN == 0:
                label = "3drac"
            specs.append(BasisSpec(label, 3, kN, 0))
        # plus anomalistic-sideband partners of the k=3 carrier
        specs.append(BasisSpec("3drac+ano",  3,  0, +1))
        specs.append(BasisSpec("3drac-ano",  3,  0, -1))
        specs.append(BasisSpec("3drac-N+ano", 3, -1, +1))
        specs.append(BasisSpec("3drac-N-ano", 3, -1, -1))
        # The 9.122 d peak is 2drac-N+ano (= drac+sid+ano); already in
        # the cross-terms block as drac-N+ano? No -- that's k=(1,-1,+1)
        # at 13.7188 d. We need the (2,-1,+1) lattice point explicitly:
        specs.append(BasisSpec("2drac-N+ano", 2, -1, +1))   # 9.1207 d
        specs.append(BasisSpec("2drac+N+ano", 2, +1, +1))   # 9.0963 d

    if use_third_harmonic:
        # Fourth-harmonic band (~6.8 d). Comes from squaring the k=2
        # lattice. Real data shows the strongest peak at 4drac-2N
        # (= 2*(2drac-N), 6.817 d, the harmonic of 13.633 d).
        for kN in (-3, -2, -1, 0, +1, +2):
            label = f"4drac{kN:+d}N".replace("+0N", "")
            if kN == 0:
                label = "4drac"
            specs.append(BasisSpec(label, 4, kN, 0))

    if use_envelopes:
        # slow envelopes (the 1+m cos(N) and 1+e cos(ano) factors
        # contribute envelope-period content too if there is any DC
        # drift in the carrier amplitude; harmless if not used)
        specs.append(BasisSpec("nodal_env",       0, +1,  0))
        specs.append(BasisSpec("anomalistic_env", 0,  0, +1))

    # De-duplicate by (k_d, k_N, k_a) -- the basis is a lattice and the
    # block-organized appends above can list the same lattice point
    # under two physical names (e.g. drac+ano appears in both the first-
    # order anomalistic block and the second-harmonic mixing block).
    seen = {}
    deduped: List[BasisSpec] = []
    for s in specs:
        key = (s.k_d, s.k_N, s.k_a)
        if key not in seen:
            seen[key] = s.label
            deduped.append(s)
    return deduped


# ----------------------------------------------------------------------------
# Linear basis matrix
# ----------------------------------------------------------------------------
def build_design_matrix(t: np.ndarray,
                        specs: List[BasisSpec],
                        detrend_degree: int = 2
                        ) -> Tuple[np.ndarray, List[str]]:
    """
    Build the OLS design matrix:
      [cos(w_1 t), sin(w_1 t),
       cos(w_2 t), sin(w_2 t), ...,
       1, t/T, (t/T)^2, ...]
    with t/T normalized to ~[-1,1] so the polynomial columns are
    well-conditioned.
    """
    t = np.asarray(t, dtype=float)
    cols = []
    names = []
    for s in specs:
        w = s.angular_freq()
        if w == 0:
            # DC (shouldn't happen here -- envelopes are nonzero)
            continue
        cols.append(np.cos(w * t))
        names.append(f"{s.label}_c")
        cols.append(np.sin(w * t))
        names.append(f"{s.label}_s")

    # Polynomial drift, normalized
    t_min, t_max = float(t.min()), float(t.max())
    span = max(t_max - t_min, 1.0)
    t_norm = 2.0 * (t - t_min) / span - 1.0   # in [-1, 1]
    for k in range(detrend_degree + 1):
        if k == 0:
            cols.append(np.ones_like(t))
            names.append("poly_0")
        else:
            cols.append(t_norm ** k)
            names.append(f"poly_{k}")

    A = np.column_stack(cols)
    return A, names


def solve_ols(A: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, float, float]:
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return coef, ss_res, r2


# ----------------------------------------------------------------------------
# Geometric back-out
# ----------------------------------------------------------------------------
def amp_phase(c_coef: float, s_coef: float) -> Tuple[float, float]:
    """
    A * cos(w t + phi)  expansion gives  A cos(phi)*cos(w t) - A sin(phi)*sin(w t).
    So if column model is c*cos(w t) + s*sin(w t), then
        A   = sqrt(c^2 + s^2)
        phi = atan2(-s, c)
    Equivalent sin form:
        A * sin(w t + psi) with psi = phi + pi/2.
    Returns (amp, phase_cos)  where phase_cos is phi above.
    """
    A = float(np.hypot(c_coef, s_coef))
    phi = float(np.arctan2(-s_coef, c_coef))
    return A, phi


@dataclass
class GeoBackout:
    """Geometric parameters reconstructed from linear amplitudes."""
    a_carrier:        float = 0.0
    d_drac:           float = 0.0     # carrier phase (cos convention)

    # nodal sideband amplitudes
    A_dN_plus:        float = 0.0     # at drac+N (period 27.1036)
    A_dN_minus:       float = 0.0     # at drac-N (period 27.3216, sidereal)
    phi_dN_plus:      float = 0.0
    phi_dN_minus:     float = 0.0
    m_N:              float = 0.0     # = (A+ + A-) / a  (locked-AM estimate)
    m_N_asymmetry:    float = 0.0     # |A+ - A-| / (A+ + A-)

    # anomalistic sideband amplitudes
    A_da_plus:        float = 0.0
    A_da_minus:       float = 0.0
    phi_da_plus:      float = 0.0
    phi_da_minus:     float = 0.0
    e_eff:            float = 0.0
    e_eff_asymmetry:  float = 0.0

    # Jacobi-Anger 3rd harmonic
    a_third_harmonic: float = 0.0
    i_eff_jacobi:     float = 0.0     # estimated from 3rd/1st ratio

    # raw amplitudes table (label -> (A, phi))
    amps: Dict[str, Tuple[float, float]] = field(default_factory=dict)


def back_out_geometric(coef: np.ndarray,
                       names: List[str]) -> GeoBackout:
    """
    Read off geometric quantities from the recovered linear coefficients.
    Convention: each component j contributes  c_j*cos(w_j t) + s_j*sin(w_j t)
    which equals  A_j * cos(w_j t + phi_j)  with A,phi from amp_phase().
    """
    # Index lookup
    idx: Dict[str, int] = {n: i for i, n in enumerate(names)}

    def amp_of(label: str) -> Tuple[float, float]:
        c = idx.get(label + "_c"); s = idx.get(label + "_s")
        if c is None or s is None:
            return (0.0, 0.0)
        return amp_phase(coef[c], coef[s])

    g = GeoBackout()
    g.amps = {}
    # carrier
    A_d, phi_d = amp_of("drac");          g.amps["drac"] = (A_d, phi_d)
    g.a_carrier = A_d
    g.d_drac = phi_d

    # nodal sidebands
    Ap, php = amp_of("drac+N");           g.amps["drac+N"] = (Ap, php)
    Am, phm = amp_of("drac-N");           g.amps["drac-N"] = (Am, phm)
    g.A_dN_plus = Ap; g.A_dN_minus = Am
    g.phi_dN_plus = php; g.phi_dN_minus = phm
    if A_d > 0:
        # AM model gives each sideband amplitude = A_d * m_N / 2
        g.m_N = (Ap + Am) / A_d
        if (Ap + Am) > 0:
            g.m_N_asymmetry = abs(Ap - Am) / (Ap + Am)

    # anomalistic sidebands
    Bp, qhp = amp_of("drac+ano");         g.amps["drac+ano"] = (Bp, qhp)
    Bm, qhm = amp_of("drac-ano");         g.amps["drac-ano"] = (Bm, qhm)
    g.A_da_plus = Bp; g.A_da_minus = Bm
    g.phi_da_plus = qhp; g.phi_da_minus = qhm
    if A_d > 0:
        g.e_eff = (Bp + Bm) / A_d
        if (Bp + Bm) > 0:
            g.e_eff_asymmetry = abs(Bp - Bm) / (Bp + Bm)

    # Jacobi-Anger 3rd harmonic
    A3, ph3 = amp_of("3drac");            g.amps["3drac"] = (A3, ph3)
    g.a_third_harmonic = A3
    if A_d > 0 and A3 > 0:
        # sin(i sin x) = 2 J_1(i) sin x + 2 J_3(i) sin 3x + ...
        # ratio J_3(i)/J_1(i)  -> for small i, ~ i^2/24
        # so i_eff ~ sqrt(24 * A3 / A_d).
        g.i_eff_jacobi = float(np.sqrt(max(0.0, 24.0 * A3 / A_d)))

    # cross terms and harmonics for completeness
    for label in ("drac+N+ano", "drac+N-ano", "drac-N+ano", "drac-N-ano",
                  "3drac+N", "3drac-N", "3drac+ano", "3drac-ano",
                  "nodal_env", "anomalistic_env"):
        A, ph = amp_of(label)
        if A > 0:
            g.amps[label] = (A, ph)

    return g


# ----------------------------------------------------------------------------
# Top-level fit
# ----------------------------------------------------------------------------
def fit_geometric_linear(y: np.ndarray, t_days: np.ndarray,
                         detrend_degree: int = 2,
                         use_third_harmonic: bool = True,
                         use_cross_terms: bool = True,
                         use_envelopes: bool = True,
                         use_perigean: bool = True,
                         use_second_harmonic: bool = True,
                         verbose: bool = True
                         ) -> Tuple[GeoBackout, np.ndarray, dict]:
    """
    Linear-basis fit. Returns (geometric_backout, model_prediction, diag).
    """
    n = len(y)
    if n < 32:
        raise ValueError(f"need at least 32 samples (got {n})")

    specs = build_basis_spec(use_third_harmonic, use_cross_terms,
                             use_envelopes, use_perigean,
                             use_second_harmonic)
    A, names = build_design_matrix(t_days, specs, detrend_degree)

    if verbose:
        print(f"[basis]  {len(specs)} sinusoidal components "
              f"({2*len(specs)} sin/cos amplitudes) + "
              f"{detrend_degree+1} polynomial = {A.shape[1]} columns")
        for s in specs:
            print(f"           {s.label:18s}  P={s.period():9.4f} d  "
                  f"(k_d={s.k_d:+d}, k_N={s.k_N:+d}, k_a={s.k_a:+d})")

    coef, ss_res, r2 = solve_ols(A, y)
    yhat = A @ coef

    if verbose:
        print(f"[ols ]  R^2 = {r2:.6f}, "
              f"RMSE = {np.sqrt(ss_res/n):.4g}")

    g = back_out_geometric(coef, names)

    # Per-band residual diagnostic
    resid = y - yhat
    band_var = _band_variance(t_days, resid)

    # Pack named coefficients for inspection/serialization
    coef_dict = {names[i]: float(coef[i]) for i in range(len(names))}

    diag = dict(
        r2 = r2,
        rmse = float(np.sqrt(ss_res/n)),
        mae = float(np.mean(np.abs(resid))),
        ss_res = ss_res,
        ss_tot = float(np.sum((y - y.mean())**2)),
        n_samples = n,
        n_basis_columns = A.shape[1],
        detrend_degree = detrend_degree,
        residual_band_variance = band_var,
        coefficients = coef_dict,
        column_names = names,
    )
    return g, yhat, diag


# ----------------------------------------------------------------------------
# Per-band residual variance diagnostic
# ----------------------------------------------------------------------------
def _band_variance(t: np.ndarray, x: np.ndarray) -> dict:
    """
    Decompose residual variance by fitting independent sin/cos pairs at
    each diagnostic period. After a successful linearized fit, ALL of
    these should be near zero except the genuine noise floor.
    """
    bands = {
        "drac":         P_DRACONIC,
        "sidereal":     P_SIDEREAL,
        "sid_mirror":   P_SID_MIRROR,
        "anomalistic":  P_ANOMALISTIC,
        "nodal":        P_NODAL,
    }
    if len(x) < 16:
        return {k: 0.0 for k in bands} | {"after_all": float(np.var(x)),
                                          "total": float(np.var(x))}

    cols = []
    for P in bands.values():
        omega = 2 * np.pi / P
        cols += [np.cos(omega * t), np.sin(omega * t)]
    A = np.column_stack(cols + [np.ones_like(t)])
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    pred_all = A @ coef
    res_all  = x - pred_all

    out = {"total": float(np.var(x)),
           "after_all": float(np.var(res_all))}
    for i, key in enumerate(bands):
        a, b = float(coef[2*i]), float(coef[2*i + 1])
        out[key] = 0.5 * (a*a + b*b)
    return out


# ----------------------------------------------------------------------------
# Synthetic test signal
# ----------------------------------------------------------------------------
def synthetic_signal(n_days: int = 8000, noise: float = 0.005,
                     seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """
    Test signal with the same generative form as the geometric model:
    a draconic carrier modulated by nodal AM + apsidal AM, plus a small
    additive sidereal-only piece to verify the residual diagnostics.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n_days, dtype=float)
    th_drac = 2 * np.pi * t / P_DRACONIC
    th_ano  = 2 * np.pi * t / P_ANOMALISTIC
    th_N    = 2 * np.pi * t / P_NODAL

    nodal_swing = np.sin(np.deg2rad(11.0) * np.sin(th_drac - 2.27))
    apsidal     = 1.0 + 0.06 * np.cos(th_ano - 2.09)
    nodal_mod   = 1.0 + 0.30 * np.cos(th_N - 1.0)
    y = 0.10 * nodal_swing * apsidal * nodal_mod

    # Small additive sidereal-only contaminant (independent of nodal AM).
    th_sid = 2 * np.pi * t / P_SIDEREAL
    y += 0.003 * np.sin(th_sid)

    y += noise * rng.standard_normal(n_days)
    return t, y


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def read_input(path: str, time_unit: str
               ) -> Tuple[np.ndarray, np.ndarray, str]:
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        y = arr
        t = np.arange(len(y), dtype=float)
        unit_used = "days (synthesized index)"
    else:
        t_raw = arr[:, 0]
        y     = arr[:, 1]
        if   time_unit == "days":    t = t_raw
        elif time_unit == "years":   t = t_raw * DAYS_PER_YEAR
        elif time_unit == "seconds": t = t_raw / 86400.0
        else: raise ValueError(f"unknown time-unit {time_unit!r}")
        unit_used = time_unit
    return t, y, unit_used


def save_outputs(out_dir: str, t: np.ndarray, y: np.ndarray,
                 yhat: np.ndarray, g: GeoBackout, diag: dict,
                 time_unit: str):
    os.makedirs(out_dir, exist_ok=True)

    bv = diag.get("residual_band_variance", {})
    var_total = float(np.var(y))
    var_resid = float(np.var(y - yhat))

    # Order amplitudes by absolute amplitude for readability
    amps_sorted = sorted(g.amps.items(),
                         key=lambda kv: -kv[1][0])

    lines = [
        "Linear-basis nodal-modulated lunisolar fit",
        "=" * 60,
        f"  n samples        : {diag['n_samples']}",
        f"  basis columns    : {diag['n_basis_columns']}",
        f"  detrend degree   : {diag['detrend_degree']}",
        f"  R^2              : {diag['r2']:.6f}",
        f"  RMSE             : {diag['rmse']:.6g}",
        f"  MAE              : {diag['mae']:.6g}",
        "",
        "Group-symmetry framing",
        "-" * 60,
        "  draconic is the k=0 carrier; sidereal (27.3216 d) and",
        "  27.1036 d emerge as the +/- nodal sidebands of draconic.",
        "  Each basis component has an independent linear amplitude --",
        "  if AM/locked-sideband holds, the +/- sidebands come out",
        "  symmetric on their own; if not, the asymmetry is the data",
        "  speaking.",
        "",
        "Sideband decomposition (sorted by amplitude)",
        "-" * 60,
        f"  {'component':<22s}  {'period (d)':>11s}  {'amplitude':>11s}  "
        f"{'phase (deg)':>11s}",
    ]
    period_of = {s.label: s.period()
                 for s in build_basis_spec(True, True, True, True, True)}
    for label, (A, phi) in amps_sorted:
        P = period_of.get(label, float("nan"))
        if P > 1e6:
            P_str = f"{P:>11.4g}"
        else:
            P_str = f"{P:>11.4f}"
        lines.append(f"  {label:<22s}  {P_str}  {A:>11.4g}  "
                     f"{np.rad2deg(phi):>11.2f}")
    lines.append("")

    lines += [
        "Geometric back-out",
        "-" * 60,
        f"  carrier amplitude a   = {g.a_carrier:.6g}",
        f"  carrier phase d_drac  = {g.d_drac:.4f} rad  "
        f"({np.rad2deg(g.d_drac):.2f} deg)",
        "",
        "  Nodal modulation (drac +/- N sidebands)",
        f"    A(drac+N, 27.1036 d) = {g.A_dN_plus:.4g}",
        f"    A(drac-N, 27.3216 d) = {g.A_dN_minus:.4g}",
        f"    m_N (locked-AM est.) = (A+ + A-)/a = {g.m_N:.4f}",
        f"    asymmetry |A+-A-|/(A++A-) = "
        f"{g.m_N_asymmetry:.4f}",
        "      (0 = perfect AM/locked sidebands; >0 = FM or "
        "independent forcing)",
        "",
        "  Apsidal modulation (drac +/- ano sidebands)",
        f"    A(drac+ano)          = {g.A_da_plus:.4g}",
        f"    A(drac-ano)          = {g.A_da_minus:.4g}",
        f"    e_eff (locked est.)  = (A+ + A-)/a = {g.e_eff:.4f}",
        f"    asymmetry            = "
        f"{g.e_eff_asymmetry:.4f}",
        "",
        "  Jacobi-Anger 3rd-harmonic check",
        f"    A(3drac)             = {g.a_third_harmonic:.4g}",
        f"    i_eff from J3/J1     = {g.i_eff_jacobi:.4f} rad  "
        f"({np.rad2deg(g.i_eff_jacobi):.3f} deg)",
        "",
        "Residual band-variance check",
        "-" * 60,
        f"  total Var(y)              : {var_total:.4g}",
        f"  Var(residual)             : {var_resid:.4g}    "
        f"({100*var_resid/var_total:.2f}% of total)",
        "",
        "  Per-period structure remaining in residual:",
        f"    draconic     27.2122 d  : {bv.get('drac',0):.4g}",
        f"    sidereal     27.3216 d  : {bv.get('sidereal',0):.4g}",
        f"    sid_mirror   27.1036 d  : {bv.get('sid_mirror',0):.4g}",
        f"    anomalistic  27.5545 d  : {bv.get('anomalistic',0):.4g}",
        f"    nodal        6798.4 d   : {bv.get('nodal',0):.4g}",
        f"    after subtracting all 5 : {bv.get('after_all',0):.4g}",
        "",
        "  These should ALL be near the noise floor after a successful",
        "  linear-basis fit. Remaining power at sidereal or sid_mirror",
        "  indicates either independent non-lunar forcing at that period",
        "  or a missing higher-order cross term.",
        "",
        "Physical-reference checks",
        "-" * 60,
        "  Lunar inclination ~ 5.14 deg",
        "  Earth obliquity   ~ 23.44 deg",
        "  Lunar e           ~ 0.0549",
        "",
        "Locked periods (days)",
        f"  P_DRAC      = {P_DRACONIC}",
        f"  P_ANO       = {P_ANOMALISTIC}",
        f"  P_NODAL     = {P_NODAL}",
        f"  P_SID  (=drac - N)  = {P_SIDEREAL}        (emergent, not free)",
        f"  P_SID' (=drac + N)  = {P_SID_MIRROR:.4f}  (emergent, not free)",
        "",
        f"  input time unit : {time_unit}",
    ]

    metrics_path = os.path.join(out_dir, "metrics.txt")
    with open(metrics_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[save]  {metrics_path}")
    print("\n".join(lines))

    # params.json
    json_path = os.path.join(out_dir, "params.json")
    out_json = {
        "geometric": asdict(g),
        "diagnostics": diag,
        "periods": dict(P_DRACONIC=P_DRACONIC,
                        P_ANOMALISTIC=P_ANOMALISTIC,
                        P_NODAL=P_NODAL,
                        P_SIDEREAL_EMERGENT=P_SIDEREAL,
                        P_SID_MIRROR_EMERGENT=P_SID_MIRROR),
    }
    # GeoBackout.amps is a Dict[str, Tuple]; tuples are JSON-fine.
    with open(json_path, "w") as f:
        json.dump(out_json, f, indent=2, default=lambda o: list(o))
    print(f"[save]  {json_path}")

    # residual.dat
    res_path = os.path.join(out_dir, "residual.dat")
    np.savetxt(res_path,
               np.column_stack([t, y, yhat, y - yhat]),
               header="t_days  y  model  residual", comments="# ")
    print(f"[save]  {res_path}")

    # fit.png
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[warn]  matplotlib unavailable; skipping fit.png")
        return

    n_show = min(len(t), 1500)
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=False)

    ax = axes[0]
    ax.plot(t[:n_show], y[:n_show], lw=0.8, color="0.3", label="data")
    ax.plot(t[:n_show], yhat[:n_show], lw=1.0, color="C3",
            label="linear-basis model")
    ax.set_xlabel("t (days)")
    ax.set_ylabel("y")
    ax.set_title(f"Linear-basis fit  -  R^2 = {diag['r2']:.4f}, "
                 f"a = {g.a_carrier:.3g}, m_N = {g.m_N:.3f}, "
                 f"e_eff = {g.e_eff:.3f}, "
                 f"asym_N = {g.m_N_asymmetry:.2f}")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(t[:n_show], (y - yhat)[:n_show], lw=0.7, color="C0")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("t (days)")
    ax.set_ylabel("residual")
    ax.set_title("residual = data - model")
    ax.grid(alpha=0.3)

    # Residual periodogram zoomed onto the 27-day band.
    # Use Hann-windowed Lomb-Scargle to suppress leakage.
    ax = axes[2]
    resid = y - yhat
    if len(resid) > 16 and _HAVE_SCIPY:
        try:
            t_f = t.astype(float)
            n = len(t_f)
            window = 0.5 - 0.5 * np.cos(2*np.pi * np.arange(n) / max(n-1, 1))
            x = (resid - resid.mean()) * window
            freqs = np.linspace(1.0/27.0, 1.0/27.1, 4000)
            angfreqs = 2 * np.pi * freqs
            pgram = lombscargle(t_f, x.astype(float),
                                angfreqs, normalize=True)
            periods = 1.0 / freqs
            ax.plot(periods, pgram, lw=0.8)
            # Perigean-band (8.85 yr) periods
            P_drac_minus_per = 1.0 / (2/P_DRACONIC - 1/P_NODAL - 1/P_ANOMALISTIC)
            P_ano_plus_N     = 1.0 / (1/P_ANOMALISTIC + 1/P_NODAL)
            P_ano_minus_N    = 1.0 / (1/P_ANOMALISTIC - 1/P_NODAL)
            markers = [
                (P_DRACONIC,        "drac",          "C3"),
                (P_SIDEREAL,        "sid (drac-N)",  "C2"),
                (P_SID_MIRROR,      "drac+N",        "C2"),
                (P_ANOMALISTIC,     "ano",           "C1"),
                (P_drac_minus_per,  "drac-per",      "C4"),
                (P_ano_plus_N,      "ano+N",         "C4"),
                (P_ano_minus_N,     "ano-N",         "C4"),
            ]
            ymax = ax.get_ylim()[1]
            for P, name, c in markers:
                if periods.min() <= P <= periods.max():
                    ax.axvline(P, color=c, lw=0.6, ls="--", alpha=0.7)
                    ax.text(P, ymax*0.95, name,
                            rotation=90, va="top", ha="right",
                            fontsize=8, color=c)
            ax.set_xlabel("period (days)")
            ax.set_ylabel("normalized power")
            ax.set_title("residual spectrum  -  all 27-d band peaks should "
                         "be near noise floor")
            ax.grid(alpha=0.3)
        except Exception as exc:
            ax.text(0.5, 0.5, f"LS error: {exc}",
                    ha="center", va="center", transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5,
                "scipy.signal.lombscargle unavailable" if not _HAVE_SCIPY
                else "not enough samples",
                ha="center", va="center", transform=ax.transAxes)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, "fit.png")
    fig.savefig(fig_path, dpi=110)
    plt.close(fig)
    print(f"[save]  {fig_path}")


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Linear-basis lunisolar geometric fitter.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("input", nargs="?", default=None,
                    help="whitespace-delimited input: 't y' or single-column y. "
                         "Omit + use --synthetic for the built-in test signal.")
    ap.add_argument("--time-unit", choices=["days", "years", "seconds"],
                    default="days",
                    help="unit of the time column (default: days)")
    ap.add_argument("--out-dir", default="geometric_out",
                    help="output directory (default: geometric_out)")
    ap.add_argument("--detrend", type=int, default=2,
                    help="degree of polynomial trend to fit alongside the "
                         "sinusoids (default 2). Use 0 to disable; 3+ for "
                         "very long spans with strong drift.")
    ap.add_argument("--no-third-harmonic", action="store_true",
                    help="disable Jacobi-Anger 3*drac terms")
    ap.add_argument("--no-cross-terms", action="store_true",
                    help="disable second-order drac+/-N+/-ano cross terms")
    ap.add_argument("--no-envelopes", action="store_true",
                    help="disable slow nodal_env / anomalistic_env terms")
    ap.add_argument("--no-perigean", action="store_true",
                    help="disable perigean (8.85 yr) sidebands: drac-per, "
                         "drac+per, ano+/-N. Enable on by default since these "
                         "explain the residual peaks at 26.985, 27.44, 27.66 d.")
    ap.add_argument("--no-second-harmonic", action="store_true",
                    help="disable half-month (~13.6-13.8 d) band: "
                         "2drac, 2drac+/-N, 2ano, 2ano+/-N. These come from "
                         "squaring/mixing first-harmonic terms.")
    ap.add_argument("--synthetic", action="store_true",
                    help="ignore --input; fit the built-in synthetic signal")
    ap.add_argument("--n-days", type=int, default=8000,
                    help="length of synthetic signal (default 8000, "
                         "covers >1 nodal cycle)")
    ap.add_argument("--noise", type=float, default=0.005,
                    help="noise sigma for synthetic signal")
    ap.add_argument("--seed",  type=int, default=42)
    args = ap.parse_args()

    if args.synthetic or args.input is None:
        print(f"[data]  using synthetic signal (n_days={args.n_days}, "
              f"noise={args.noise})")
        t, y = synthetic_signal(n_days=args.n_days,
                                noise=args.noise, seed=args.seed)
        time_unit = "days (synthetic)"
        os.makedirs(args.out_dir, exist_ok=True)
        np.savetxt(os.path.join(args.out_dir, "synthetic_input.dat"),
                   np.column_stack([t, y]),
                   header="t_days  y", comments="# ")
    else:
        t, y, time_unit = read_input(args.input, args.time_unit)
        print(f"[data]  read {len(y)} samples from {args.input}  "
              f"(time_unit={time_unit})")
        if len(t) > 1:
            span = t.max() - t.min()
            print(f"        span = {span:.1f} days "
                  f"= {span/DAYS_PER_YEAR:.2f} years "
                  f"({span/P_NODAL:.2f} nodal cycles)")

    print(f"[fit ]  linear-basis OLS, detrend degree = {args.detrend}")
    g, yhat, diag = fit_geometric_linear(
        y, t,
        detrend_degree=args.detrend,
        use_third_harmonic=not args.no_third_harmonic,
        use_cross_terms=not args.no_cross_terms,
        use_envelopes=not args.no_envelopes,
        use_perigean=not args.no_perigean,
        use_second_harmonic=not args.no_second_harmonic,
        verbose=True,
    )
    save_outputs(args.out_dir, t, y, yhat, g, diag, time_unit)


if __name__ == "__main__":
    main()
