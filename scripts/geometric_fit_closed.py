#!/usr/bin/env python3
"""
geometric_fit_closed.py
=======================

Constrained nonlinear closed-form fitter for the lunisolar signal.
Companion to geometric_fit.py (which solves the unconstrained linear
expansion on the (k_d, k_N, k_a) lattice).

The original closed form
------------------------

    y(t) = a * sin( i_eff sin(theta_d + d_d) )
             * (1 + e_eff cos(theta_a + d_a))
             * (1 + m_N   cos(theta_N + d_N))
           + b0

has only 8 parameters and predicts strict relationships between
sideband amplitudes:

    A(drac +/- N)        = a * J_1(i) * m_N / 2          (locked equal)
    A(drac +/- ano)      = a * J_1(i) * e_eff / 2        (locked equal)
    A(drac +/- N +/- ano)= a * J_1(i) * m_N * e_eff / 4
    A(3drac +/- ...)     = a * J_3(i) ...
    A(2drac, 2ano, ...)  = 0                              (forbidden!)

The linear-basis fit revealed that the *forbidden* even-harmonic peaks
(13.633, 13.606, 13.777, 6.817, ... d) are present in real data with
substantial amplitude. So the strict closed form above CANNOT match
the data: it produces only odd k_d harmonics by Jacobi-Anger.

Two closed forms supported here
-------------------------------

1. ``pure`` -- the 11-parameter physical-geometry form:

        y(t) = a * [ c_0 + sin(i sin(th_d + d_drac)) ]
                 * (1 + e_eff cos(th_a + d_ano))
                 * (1 + m_N   cos(th_N + d_N))
             + a_sid * cos(th_sid + d_sid)
                     * (1 + e_eff cos(th_a + d_ano))
             + b0

   The leading c_0 lets a *standalone* anomalistic peak emerge at
   P=27.55 d as a*c_0*e_eff*cos(th_a+d_ano), which the original AC-only
   carrier `sin(i*sin th_d)` could never produce. The sidereal block
   (a_sid, d_sid) is an additive sidereal harmonic at P=27.3216 d (sid =
   drac - N) modulated by the same anomalistic envelope, so a SINGLE
   amp+phase pair simultaneously emits both the sidereal peak and the
   anomalistic+sidereal cross-term at P=13.719 d. Both are physically
   essential because m_N tends to collapse to ~0 in real data, which
   kills the sideband structure that the carrier alone would otherwise
   create. (The retired 3*sid term, which assumed the strong 'sid_mirror'
   at drac+N, was a worse model -- the +N mirror is weak, the -N
   sidereal is strong.)

2. ``extended`` -- pure + an explicit even-harmonic generator (i2,
   d_drac2) inside the inclination factor:

        swing -> c_0 + sin(i*sin(th_d+d_d)) + i2 * sin(2*(th_d+d_d) + d_d2)

   13 parameters. The 'i2' term breaks the Jacobi-Anger no-2*drac
   restriction and lets a free 13.6 d band emerge. (Compound tides
   typically explain that band better -- see --compound-tides lunar.)

Compound tides (--compound-tides lunar / --strict-orbital)
----------------------------------------------------------

The long-period tidal spectrum has TWO independent origins that the
orbital-product model conflates:

  1. ORBITAL band (~27 d): direct gravitational geometry of the
     draconic, anomalistic, and nodal cycles -- the strict closed
     form is the right structure here.

  2. COMPOUND-TIDE band (~6-15 d): NOT geometric harmonics, but
     non-linear beats between fast (sub-daily) constituents (K1, O1,
     M2, S2, ...) created by friction in shallow seas and leaked back
     into the open ocean. They land at:

        Mf   @ 13.661 d   (K1 - O1, lunar fortnightly)
        Mt   @  9.133 d   (3 M2 - 2 S2-ish, ter-mensal lunar)
        MSqm @  6.831 d   (compound quarter-monthly)

  This means the 13.6 d and 6.8 d residual peaks the linear basis
  identifies as '2*drac' / '4*drac-2N' may be misattributed: they
  could be Mf and MSqm, sharing the band by spectral coincidence.

When --compound-tides lunar is enabled, an ADDITIVE compound-tide
module is appended to the orbital base:

        y(t) = orbital_base(t) + compound(t) + b0
        compound(t) = A_Mf  cos(2*pi t/13.661 + phi_Mf)
                    + A_Mt  cos(2*pi t/9.133  + phi_Mt)
                    + A_MSqm cos(2*pi t/6.831 + phi_MSqm)

Periods are LOCKED to the canonical compound-tide values; only
amplitude and phase are free (6 extra params total).

Combine with --strict-orbital to drop the i2 (even-harmonic) term
from the extended orbital model. This is the clean test of the
question 'is the 13.6 d band a geometric Bessel artifact (i2 needed)
or compound-tide contamination (compound module needed, i2 redundant)?'
If the strict-orbital + compound combination beats extended without
compound, then the i2 term in the extended model was dressing up
compound-tide signal as if it were orbital geometry.

Flag matrix (all combine with --mix-semiannual / --mix-evection):
  --model pure                                  : 11 params, baseline
  --model extended                              : 13 params, baseline
  --model extended  --compound-tides lunar      : 19 params, both modules
  --model extended  --compound-tides lunar  --strict-orbital
                                                : 17 params, clean
                                                  separation: orbital is
                                                  pure 11-param form +
                                                  compound is 6 free.
  --model pure      --compound-tides lunar      : 17 params, equivalent
                                                  to the above (pure is
                                                  already strict)

Mixing carriers (--mix-semiannual / --mix-evection)
---------------------------------------------------

When the upstream signal is generated by stroboscopic sampling at the
annual rate (or a harmonic of it), low-frequency carriers like the
semi-annual solar tide (~182.6 d) and the classical lunar evection
(~31.81 d) are folded down to DC and re-emerge as ALIASES sitting next
to the lunar fundamentals. Symptoms in the residual:

    9.083 d  <-  Mst tidal factor at 9.557 d, aliased by ~1/(half-year)
    27.093 d <-  classical evection at 31.81 d, aliased by ~1/year
    asymmetric drac-N vs drac+N sidebands -- standalone sid/ano peaks

If those carriers are real but unobservable post-strobe, the optimizer
cannot simply fit them as additive components: you have to *include
them in the forward mixing* (multiplicative) so the lunar amplitudes
are no longer biased to absorb the alias spillage, and *remove them
before scoring* so R^2 is judged against the de-aliased lunar skeleton.

Forward model with mixing carriers active:

    M(t)  = 1 + alpha_sa * cos(2*pi*t/T_sa + phi_sa)        # semi-annual
                + alpha_ev * cos(2*pi*t/T_ev + phi_ev)      # evection-aliased
    y(t)  = base(t) * M(t) + b0

Residual (subtraction filter):

    de_aliased_data  = y - base(t) * (M(t) - 1)
    residual         = de_aliased_data - (base(t) + b0)
                     = y - base(t) * M(t) - b0

That is, residual is computed against the FULL forward model -- which
is the correct thing to score (the optimizer must explain y, not the
carrier-stripped version of itself). The point of the subtraction is
that the physics interpretation of the lunar amplitudes is now clean:
they describe the unmodulated lunar geometry, while the alpha_sa,
alpha_ev numbers carry the aliasing-carrier strength explicitly.

Flags:
    --mix-semiannual  : enables alpha_sa/T_sa/phi_sa  (3 extra params)
    --mix-evection    : enables alpha_ev/T_ev/phi_ev  (3 extra params)

Both can be combined with --model {pure,extended}. T_sa is bounded
[170, 195] d and T_ev [29, 35] d so the optimizer cannot escape into
unrelated low-frequency content.

Bootstrap
---------
Both variants are initialized from the linear-basis fit's recovered
amplitudes -- specifically:

    a      <- amplitude at the (1,0,0) drac column
    d_drac <- phase at (1,0,0)
    m_N    <- 2 * mean(A(drac+N), A(drac-N)) / a
    e_eff  <- 2 * mean(A(drac+ano), A(drac-ano)) / a
    d_N, d_ano  <- inferred from sideband phase combinations
    i_eff  <- from the J_3/J_1 ratio of (3drac)/(drac) amplitudes,
              when present, else seeded at 0.1 rad.

This collapses the multi-modal nonlinear landscape that plagued the
old multistart fitter to a single warm-started TRF refinement. It
runs in <1 s on tens of thousands of samples.

Usage
-----
    python geometric_fit_closed.py daily.dat --time-unit years
    python geometric_fit_closed.py --synthetic
    python geometric_fit_closed.py daily.dat --model pure
    python geometric_fit_closed.py daily.dat --model extended  # default
    python geometric_fit_closed.py daily.dat --mix-semiannual
    python geometric_fit_closed.py daily.dat --mix-semiannual --mix-evection

The script writes the same set of outputs as geometric_fit.py
(metrics.txt, params.json, residual.dat, fit.png).
"""

from __future__ import annotations
import argparse
import json
import os
from dataclasses import dataclass, asdict
from typing import Tuple

import numpy as np

try:
    from scipy.optimize import least_squares
except ImportError:
    raise SystemExit(
        "scipy is required: pip install scipy\n"
        "(Same environment as lunisolar_pysr.py / geometric_fit.py.)"
    )

# Reuse infrastructure from the linear-basis fitter
import geometric_fit as gf


# ============================================================================
# Closed-form models
# ============================================================================
# Default canonical periods for the mixing carriers, in days.
T_SEMIANNUAL_DEFAULT = 182.6211   # tropical half-year
T_EVECTION_DEFAULT   = 31.811     # classical lunar evection

# Canonical compound-tide periods (days) -- LOCKED in the fitter,
# only amplitude and phase are free.
T_Mf   = 13.661   # K1 - O1, lunar fortnightly
T_Mt   = 9.133    # ter-mensal lunar
T_MSqm = 6.831    # compound quarter-monthly


@dataclass
class MixSpec:
    """Which mixing carriers are active."""
    semiannual: bool = False
    evection:   bool = False

    @property
    def n_extra(self) -> int:
        return 3 * int(self.semiannual) + 3 * int(self.evection)

    @property
    def label(self) -> str:
        bits = []
        if self.semiannual: bits.append("semiannual")
        if self.evection:   bits.append("evection")
        return "+".join(bits) if bits else "none"


@dataclass
class CompoundSpec:
    """Which compound-tide additive terms are active."""
    lunar: bool = False           # adds Mf + Mt + MSqm (6 params total)

    @property
    def n_extra(self) -> int:
        return 6 * int(self.lunar)

    @property
    def label(self) -> str:
        return "lunar" if self.lunar else "none"

    @property
    def names(self) -> Tuple[str, ...]:
        return ("Mf", "Mt", "MSqm") if self.lunar else ()


@dataclass
class StrictSpec:
    """Whether to suppress extended-model orbital additions."""
    orbital: bool = False         # drop i2 in extended model


@dataclass
class HarmonicSpec:
    """Explicit free harmonic amplitudes on the lattice.

    The closed-form product expansion produces specific amplitude *ratios*
    between sidebands (e.g. drac+N == drac-N forced to be equal). On real
    data these ratios are violated. This module adds explicit free
    amplitudes for the lattice points where the gap is largest, breaking
    the strict ratio constraints while keeping the model interpretable.

    Levels:
      'none' -- no extra free harmonics
      'k3'   -- adds 3*drac block: (3,0,0), (3,-1,0), (3,-2,0), (3,-3,0)
                plus the 2drac block (2,0,0), (2,-1,0)
                plus 2ano (0,0,2) -- 7 amp/phase pairs = 14 params
      'full' -- 'k3' + asymmetric carrier sidebands (drac+N, drac-N
                additive corrections allowing amplitude asymmetry)
                plus drac-N+ano, 2drac+N+ano, ano-N
                -- 12 amp/phase pairs = 24 params
    """
    level: str = "none"           # 'none' | 'k3' | 'full'

    @property
    def n_extra(self) -> int:
        if self.level == "k3":
            return 14
        if self.level == "full":
            return 24
        return 0

    @property
    def label(self) -> str:
        return self.level


# Lattice points covered at each level (kd, kN, ka).
# Each contributes one (A, phi) free pair.
_HARM_K3_LATTICE = (
    # 3*drac block (the 9-day band)
    ("H_3drac",      3,  0,  0),    #  9.0707 d
    ("H_3drac_m1N",  3, -1,  0),    #  9.0829 d  <-- biggest gap
    ("H_3drac_m2N",  3, -2,  0),    #  9.0950 d
    ("H_3drac_m3N",  3, -3,  0),    #  9.1072 d
    # 2*drac block (the 13.6 d band)
    ("H_2drac",      2,  0,  0),    # 13.6061 d
    ("H_2drac_m1N",  2, -1,  0),    # 13.6334 d
    # 2*ano peak (the 13.78 d secondary)
    ("H_2ano",       0,  0,  2),    # 13.7773 d
)
_HARM_FULL_EXTRA_LATTICE = (
    # asymmetric carrier sidebands (additive corrections):
    ("H_drac_p1N",   1, +1,  0),    # 27.1037 d
    ("H_drac_m1N",   1, -1,  0),    # 27.3216 d
    # cross-band gaps
    ("H_drac_m1N_p1a", 1, -1, +1),   # 13.7188 d
    ("H_2drac_p1N_p1a", 2, +1, +1), #  9.0963 d
    ("H_ano_m1N",    0, -1, +1),    # 27.6666 d
)


def _harm_lattice(level: str):
    if level == "k3":
        return _HARM_K3_LATTICE
    if level == "full":
        return _HARM_K3_LATTICE + _HARM_FULL_EXTRA_LATTICE
    return ()


def _harm_omega(kd: int, kN: int, ka: int) -> float:
    """Angular frequency of lattice point (kd, kN, ka), in 1/day."""
    return 2*np.pi*(kd/gf.P_DRACONIC + kN/gf.P_NODAL + ka/gf.P_ANOMALISTIC)


@dataclass
class ClosedParams:
    """Parameters of the closed-form models.

    Pure model (v4): 11 parameters.
        y = a * [c_0 + sin(i*sin(th_d + d_drac))]
              * (1 + e_eff * cos(th_a + d_ano))
              * (1 + m_N   * cos(th_N + d_N))
          + a_sid * cos(th_sid + d_sid)
                  * (1 + e_eff * cos(th_a + d_ano))
          + b0

    The leading c_0 inside the inclination factor lets the standalone
    anomalistic peak (P=27.55 d) emerge as a*c_0*e_eff*cos(th_a+d_ano).
    The sidereal block (a_sid, d_sid) provides direct sidereal energy at
    P=27.32 d, and -- through the (1+e*cos th_a) modulation -- it also
    spawns the anomalistic-sidereal cross-term at P=13.719 d. Both are
    needed because m_N tends to collapse to ~0 in real data, which kills
    the sideband structure that the carrier alone would otherwise create.
    """
    # Pure-model 11 parameters
    a:       float = 0.0
    c_0:     float = 0.0       # DC level inside inclination factor
    i_eff:   float = 0.0       # rad
    e_eff:   float = 0.0
    m_N:     float = 0.0
    d_drac:  float = 0.0
    d_ano:   float = 0.0
    d_N:     float = 0.0
    a_sid:   float = 0.0       # sidereal block amplitude
    d_sid:   float = 0.0       # sidereal block phase
    b0:      float = 0.0
    # Extended-model additions (only i2 remains; a_3sid retired)
    i2:      float = 0.0       # even-harmonic amplitude (relative to a)
    d_drac2: float = 0.0       # 2*drac phase offset
    # Mixing-carrier additions (semi-annual)
    alpha_sa: float = 0.0      # multiplicative depth (dimensionless)
    T_sa:     float = T_SEMIANNUAL_DEFAULT
    phi_sa:   float = 0.0
    # Mixing-carrier additions (evection)
    alpha_ev: float = 0.0
    T_ev:     float = T_EVECTION_DEFAULT
    phi_ev:   float = 0.0
    # Compound-tide additions (additive; locked periods)
    A_Mf:    float = 0.0
    phi_Mf:  float = 0.0
    A_Mt:    float = 0.0
    phi_Mt:  float = 0.0
    A_MSqm:  float = 0.0
    phi_MSqm: float = 0.0
    # Free harmonic amplitudes (k3 level)
    A_H_3drac:        float = 0.0
    phi_H_3drac:      float = 0.0
    A_H_3drac_m1N:    float = 0.0
    phi_H_3drac_m1N:  float = 0.0
    A_H_3drac_m2N:    float = 0.0
    phi_H_3drac_m2N:  float = 0.0
    A_H_3drac_m3N:    float = 0.0
    phi_H_3drac_m3N:  float = 0.0
    A_H_2drac:        float = 0.0
    phi_H_2drac:      float = 0.0
    A_H_2drac_m1N:    float = 0.0
    phi_H_2drac_m1N:  float = 0.0
    A_H_2ano:         float = 0.0
    phi_H_2ano:       float = 0.0
    # Free harmonic amplitudes (full level extras)
    A_H_drac_p1N:        float = 0.0
    phi_H_drac_p1N:      float = 0.0
    A_H_drac_m1N:        float = 0.0
    phi_H_drac_m1N:      float = 0.0
    A_H_drac_m1N_p1a:    float = 0.0
    phi_H_drac_m1N_p1a:  float = 0.0
    A_H_2drac_p1N_p1a:   float = 0.0
    phi_H_2drac_p1N_p1a: float = 0.0
    A_H_ano_m1N:         float = 0.0
    phi_H_ano_m1N:       float = 0.0


# Sidereal angular speed: sid = drac - N, so omega_sid = omega_drac - omega_N.
def _base_pure(t: np.ndarray, p: ClosedParams) -> np.ndarray:
    th_d   = 2*np.pi*t/gf.P_DRACONIC
    th_a   = 2*np.pi*t/gf.P_ANOMALISTIC
    th_N   = 2*np.pi*t/gf.P_NODAL
    th_sid = th_d - th_N                                 # P=27.3216 d
    swing   = p.c_0 + np.sin(p.i_eff * np.sin(th_d + p.d_drac))
    apsidal = 1.0 + p.e_eff * np.cos(th_a + p.d_ano)
    nodal   = 1.0 + p.m_N   * np.cos(th_N + p.d_N)
    base    = p.a * swing * apsidal * nodal
    # Sidereal block (replaces a_3sid). Modulated by the same anomalistic
    # envelope so that the anomalistic+sidereal cross-term at 13.719 d
    # emerges from a single additional pair of parameters.
    sid_block = p.a_sid * np.cos(th_sid + p.d_sid) * apsidal
    return base + sid_block


def _base_extended(t: np.ndarray, p: ClosedParams,
                   strict: StrictSpec | None = None) -> np.ndarray:
    th_d   = 2*np.pi*t/gf.P_DRACONIC
    th_a   = 2*np.pi*t/gf.P_ANOMALISTIC
    th_N   = 2*np.pi*t/gf.P_NODAL
    th_sid = th_d - th_N
    th_d_phi = th_d + p.d_drac
    if strict is not None and strict.orbital:
        # Drop the i2 addition: extended collapses to pure.
        swing = p.c_0 + np.sin(p.i_eff * np.sin(th_d_phi))
    else:
        swing = (p.c_0
                 + np.sin(p.i_eff * np.sin(th_d_phi))
                 + p.i2 * np.sin(2.0 * th_d_phi + p.d_drac2))
    apsidal = 1.0 + p.e_eff * np.cos(th_a + p.d_ano)
    nodal   = 1.0 + p.m_N   * np.cos(th_N + p.d_N)
    base    = p.a * swing * apsidal * nodal
    # Sidereal block (always active; subsumes the retired a_3sid term).
    sid_block = p.a_sid * np.cos(th_sid + p.d_sid) * apsidal
    return base + sid_block


def _compound_term(t: np.ndarray, p: ClosedParams,
                   compound: CompoundSpec | None) -> np.ndarray:
    """Additive compound-tide module (locked canonical periods)."""
    out = np.zeros_like(t, dtype=float)
    if compound is None or not compound.lunar:
        return out
    out = out + p.A_Mf   * np.cos(2*np.pi*t/T_Mf   + p.phi_Mf)
    out = out + p.A_Mt   * np.cos(2*np.pi*t/T_Mt   + p.phi_Mt)
    out = out + p.A_MSqm * np.cos(2*np.pi*t/T_MSqm + p.phi_MSqm)
    return out


def _harmonic_term(t: np.ndarray, p: ClosedParams,
                   harm: HarmonicSpec | None) -> np.ndarray:
    """Additive free-harmonic block on the (kd, kN, ka) lattice."""
    out = np.zeros_like(t, dtype=float)
    if harm is None or harm.level == "none":
        return out
    for label, kd, kN, ka in _harm_lattice(harm.level):
        A = getattr(p, "A_" + label)
        phi = getattr(p, "phi_" + label)
        omega = _harm_omega(kd, kN, ka)
        out = out + A * np.cos(omega * t + phi)
    return out


def _mix_envelope(t: np.ndarray, p: ClosedParams, mix: MixSpec) -> np.ndarray:
    """Returns M(t) - 1: the multiplicative envelope minus the DC level."""
    env = np.zeros_like(t, dtype=float)
    if mix.semiannual:
        env = env + p.alpha_sa * np.cos(2*np.pi*t/p.T_sa + p.phi_sa)
    if mix.evection:
        env = env + p.alpha_ev * np.cos(2*np.pi*t/p.T_ev + p.phi_ev)
    return env


def evaluate(t: np.ndarray, p: ClosedParams, model: str,
             mix: MixSpec | None = None,
             compound: CompoundSpec | None = None,
             strict: StrictSpec | None = None,
             harm: HarmonicSpec | None = None) -> np.ndarray:
    """Forward model: y_hat = base(t) * (1 + envelope(t)) + compound + harmonics + b0."""
    if model == "pure":
        base = _base_pure(t, p)
    elif model == "extended":
        base = _base_extended(t, p, strict)
    else:
        raise ValueError(f"unknown model {model!r}")
    if mix is not None and (mix.semiannual or mix.evection):
        env = _mix_envelope(t, p, mix)
        out = base * (1.0 + env)
    else:
        out = base
    out = out + _compound_term(t, p, compound)
    out = out + _harmonic_term(t, p, harm)
    return out + p.b0


def evaluate_demixed(t: np.ndarray, p: ClosedParams, model: str,
                     mix: MixSpec | None = None,
                     compound: CompoundSpec | None = None,
                     strict: StrictSpec | None = None,
                     harm: HarmonicSpec | None = None) -> np.ndarray:
    """
    The de-aliased lunar skeleton: forward model with the mixing
    envelope subtracted (filter step). Compound and harmonic modules
    are included if active. Equivalent to evaluate() with mix disabled,
    evaluated at the mix-aware best-fit parameters.
    """
    if model == "pure":
        base = _base_pure(t, p)
    else:
        base = _base_extended(t, p, strict)
    return (base + _compound_term(t, p, compound)
            + _harmonic_term(t, p, harm) + p.b0)


def demixed_data(t: np.ndarray, y: np.ndarray, p: ClosedParams,
                 model: str, mix: MixSpec,
                 strict: StrictSpec | None = None) -> np.ndarray:
    """
    Apply the subtraction filter to the data:
        y_dealiased = y - base(t) * envelope(t)
    so the alias-carrying sidebands are stripped before residual
    interpretation. Score-relevant residual is still y - full_model.
    """
    if not (mix.semiannual or mix.evection):
        return y
    if model == "pure":
        base = _base_pure(t, p)
    else:
        base = _base_extended(t, p, strict)
    env = _mix_envelope(t, p, mix)
    return y - base * env


# ============================================================================
# Bootstrap from linear-basis fit
# ============================================================================
def bootstrap_from_linear(y: np.ndarray, t: np.ndarray,
                          verbose: bool = True
                          ) -> Tuple[ClosedParams, dict]:
    """
    Run the linear-basis fit (from geometric_fit.py) and read off
    initial estimates for the closed-form parameters.
    """
    if verbose:
        print("[boot]  running linear-basis OLS to seed closed form...")
    g_lin, yhat_lin, diag_lin = gf.fit_geometric_linear(
        y, t, detrend_degree=2,
        use_third_harmonic=True, use_cross_terms=True,
        use_envelopes=True, use_perigean=True,
        use_second_harmonic=True,
        verbose=False,
    )

    p = ClosedParams()
    p.a      = float(g_lin.a_carrier)
    # In our cos-convention, model_pure(a*sin(th+d)) at small i ~ a*i*sin(th+d).
    # The linear basis 'drac' column had phase phi (cos-form):
    # A cos(w t + phi) = A sin(w t + phi + pi/2). So d_drac for the
    # sin-form model = phi + pi/2.
    p.d_drac = float(g_lin.d_drac + np.pi/2)
    p.b0     = 0.0

    # m_N from the locked-AM estimate already computed
    p.m_N   = max(0.0, min(1.0, float(g_lin.m_N)))
    # Initial nodal phase: average of the two sideband phases (cos-conv).
    if g_lin.A_dN_plus + g_lin.A_dN_minus > 0:
        p.d_N = 0.5 * (g_lin.phi_dN_plus + g_lin.phi_dN_minus)
    else:
        p.d_N = 0.0

    p.e_eff = max(0.0, min(0.9, float(g_lin.e_eff)))
    if g_lin.A_da_plus + g_lin.A_da_minus > 0:
        p.d_ano = 0.5 * (g_lin.phi_da_plus + g_lin.phi_da_minus)
    else:
        p.d_ano = 0.0

    # i_eff from the Bessel ratio (already computed in back-out)
    p.i_eff = float(g_lin.i_eff_jacobi) if g_lin.i_eff_jacobi > 0 else 0.1
    p.i_eff = max(0.01, min(np.pi/2, p.i_eff))

    # New pure-model seeds: c_0 (DC of inclination factor) and the
    # sidereal block (a_sid, d_sid).
    #
    # c_0 has no direct linear-basis observable; seed it small but nonzero
    # so the optimizer's gradient is nonzero in that direction.
    p.c_0 = 0.05
    # a_sid, d_sid: read directly from the linear-basis amplitude at the
    # sidereal lattice point (drac-N), which lives at P=27.3216 d.
    A_sid, phi_sid = g_lin.amps.get("drac-N", (0.0, 0.0))
    p.a_sid = float(A_sid)
    p.d_sid = float(phi_sid)

    # Extended-model seeds
    # i2: from amplitude at (k_d=2, k_N=0, k_a=0) i.e. '2drac' column,
    # rescaled by carrier a.
    A_2d, phi_2d = g_lin.amps.get("2drac", (0.0, 0.0))
    p.i2      = (A_2d / p.a) if p.a > 0 else 0.0
    # in extended model, even-harmonic uses sin(2*(th_d+d_d) + d_d2);
    # cos basis at 2*drac means amplitude phase = phi_2d (cos-conv).
    # sin form -> add pi/2 and subtract 2*d_drac.
    p.d_drac2 = float(phi_2d + np.pi/2 - 2*p.d_drac)

    # Mixing-carrier seeds: small nonzero alpha so the optimizer's
    # gradient is nonzero in those directions. Period seeds at the
    # canonical values; phase 0.
    p.alpha_sa = 0.05
    p.T_sa     = T_SEMIANNUAL_DEFAULT
    p.phi_sa   = 0.0
    p.alpha_ev = 0.05
    p.T_ev     = T_EVECTION_DEFAULT
    p.phi_ev   = 0.0

    # Compound-tide seeds: read directly from the linear-basis amplitudes
    # at the corresponding lattice points, since those are what currently
    # absorb the compound-tide signal. (Mf=13.66 is near 2drac=13.61;
    # Mt=9.13 is near drac+2sid=9.094; MSqm=6.83 is near 4drac-2N=6.817.)
    # If those columns aren't present in the basis, default to 0.
    A_2d, _   = g_lin.amps.get("2drac",      (0.0, 0.0))
    A_dr2s, _ = g_lin.amps.get("drac+2sid",  (0.0, 0.0))
    A_4d2N, _ = g_lin.amps.get("4drac-2N",   (0.0, 0.0))
    p.A_Mf    = float(A_2d)
    p.phi_Mf  = 0.0
    p.A_Mt    = float(A_dr2s)
    p.phi_Mt  = 0.0
    p.A_MSqm  = float(A_4d2N)
    p.phi_MSqm = 0.0

    if verbose:
        print(f"[boot]  seeded from linear: a={p.a:.4g}, c_0={p.c_0:.4g}, "
              f"i_eff={np.rad2deg(p.i_eff):.2f} deg, "
              f"e_eff={p.e_eff:.4f}, m_N={p.m_N:.4f}, "
              f"a_sid={p.a_sid:.4g}, i2={p.i2:.4g}")

    return p, dict(
        linear_r2 = diag_lin["r2"],
        linear_rmse = diag_lin["rmse"],
        linear_amps = {k: list(v) for k, v in g_lin.amps.items()},
    )


# ============================================================================
# Bounded TRF refinement
# ============================================================================
# Core 11-vector  -- common to pure & extended, mix-independent.
# Includes c_0 (DC inside inclination factor) and the sidereal block
# (a_sid, d_sid). The sidereal block replaces the retired a_3sid term.
_CORE_NAMES = ("a", "c_0", "i_eff", "e_eff", "m_N",
               "d_drac", "d_ano", "d_N",
               "a_sid", "d_sid", "b0")
_EXT_NAMES  = ("i2", "d_drac2")
_SA_NAMES   = ("alpha_sa", "T_sa", "phi_sa")
_EV_NAMES   = ("alpha_ev", "T_ev", "phi_ev")
_CMP_NAMES  = ("A_Mf", "phi_Mf", "A_Mt", "phi_Mt", "A_MSqm", "phi_MSqm")

_CORE_LB = np.array([1e-9,  -1.0, 1e-3, 0.0, 0.0,
                     -10*np.pi, -10*np.pi, -10*np.pi,
                     -np.inf, -10*np.pi, -np.inf])
_CORE_UB = np.array([np.inf, 1.0, np.pi/2, 0.9, 1.5,
                      10*np.pi,  10*np.pi,  10*np.pi,
                      np.inf, 10*np.pi, np.inf])
_EXT_LB  = np.array([-1.0, -10*np.pi])
_EXT_UB  = np.array([ 1.0,  10*np.pi])
# Mixing carriers: alpha bounded ~[-0.9, 0.9] so 1+env stays positive,
# T bounded tightly so the optimizer can't escape.
_SA_LB = np.array([-0.9, 170.0, -10*np.pi])
_SA_UB = np.array([ 0.9, 195.0,  10*np.pi])
_EV_LB = np.array([-0.9,  29.0, -10*np.pi])
_EV_UB = np.array([ 0.9,  35.0,  10*np.pi])
# Compound-tide additive amplitudes: free amplitude (any sign), free
# phase. Periods are LOCKED (not optimized).
_CMP_LB = np.array([-np.inf, -10*np.pi]*3)
_CMP_UB = np.array([ np.inf,  10*np.pi]*3)
# Free harmonic amplitudes: free amp (positive), free phase. Periods locked.
# Names are auto-derived from the lattice tuples below.
def _harm_names_for_level(level: str) -> Tuple[str, ...]:
    out: Tuple[str, ...] = ()
    for label, _, _, _ in _harm_lattice(level):
        out = out + ("A_" + label, "phi_" + label)
    return out

def _harm_lb_ub_for_level(level: str) -> Tuple[np.ndarray, np.ndarray]:
    n = len(_harm_lattice(level))
    if n == 0:
        return np.array([]), np.array([])
    lb = np.tile(np.array([0.0, -10*np.pi]), n)
    ub = np.tile(np.array([np.inf, 10*np.pi]), n)
    return lb, ub


def _layout(model: str, mix: MixSpec,
            compound: CompoundSpec | None = None,
            strict:   StrictSpec   | None = None,
            harm:     HarmonicSpec | None = None,
            ) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray]:
    """Return (name-tuple, lower-bound, upper-bound) for the active layout."""
    if compound is None: compound = CompoundSpec()
    if strict   is None: strict   = StrictSpec()
    if harm     is None: harm     = HarmonicSpec()
    names = _CORE_NAMES
    lb, ub = _CORE_LB, _CORE_UB
    if model == "extended" and not strict.orbital:
        # Only include i2 as free param if extended AND not strict.
        names = names + _EXT_NAMES
        lb = np.concatenate([lb, _EXT_LB])
        ub = np.concatenate([ub, _EXT_UB])
    elif model not in ("pure", "extended"):
        raise ValueError(f"unknown model {model!r}")
    if mix.semiannual:
        names = names + _SA_NAMES
        lb = np.concatenate([lb, _SA_LB])
        ub = np.concatenate([ub, _SA_UB])
    if mix.evection:
        names = names + _EV_NAMES
        lb = np.concatenate([lb, _EV_LB])
        ub = np.concatenate([ub, _EV_UB])
    if compound.lunar:
        names = names + _CMP_NAMES
        lb = np.concatenate([lb, _CMP_LB])
        ub = np.concatenate([ub, _CMP_UB])
    if harm.level != "none":
        h_names = _harm_names_for_level(harm.level)
        h_lb, h_ub = _harm_lb_ub_for_level(harm.level)
        names = names + h_names
        lb = np.concatenate([lb, h_lb])
        ub = np.concatenate([ub, h_ub])
    return names, lb, ub


def _pack(p: ClosedParams, names: Tuple[str, ...]) -> np.ndarray:
    return np.array([float(getattr(p, n)) for n in names])


def _unpack(x: np.ndarray, names: Tuple[str, ...],
            template: ClosedParams | None = None) -> ClosedParams:
    p = ClosedParams() if template is None else ClosedParams(**asdict(template))
    for name, val in zip(names, x):
        setattr(p, name, float(val))
    return p


def fit_closed(y: np.ndarray, t: np.ndarray, p0: ClosedParams,
               model: str = "extended",
               mix: MixSpec | None = None,
               compound: CompoundSpec | None = None,
               strict: StrictSpec | None = None,
               harm: HarmonicSpec | None = None,
               verbose: bool = True
               ) -> Tuple[ClosedParams, dict]:
    if mix is None:      mix = MixSpec()
    if compound is None: compound = CompoundSpec()
    if strict is None:   strict = StrictSpec()
    if harm is None:     harm = HarmonicSpec()
    n = len(y)
    names, lb, ub = _layout(model, mix, compound, strict, harm)

    x0 = _pack(p0, names)
    # Snap initial guess into bounds (TRF requires strict feasibility)
    x0 = np.clip(x0, lb + 1e-12, ub - 1e-12)

    def residual(x):
        return evaluate(t, _unpack(x, names, p0), model,
                        mix, compound, strict, harm) - y

    res = least_squares(residual, x0, bounds=(lb, ub),
                        method="trf", max_nfev=8000, xtol=1e-12)
    p_fit = _unpack(res.x, names, p0)

    # Wrap phases into [0, 2pi)
    phase_names = ["d_drac", "d_ano", "d_N", "d_sid", "d_drac2",
                   "phi_sa", "phi_ev",
                   "phi_Mf", "phi_Mt", "phi_MSqm"]
    for label, _, _, _ in _harm_lattice(harm.level):
        phase_names.append("phi_" + label)
    for name in phase_names:
        v = getattr(p_fit, name)
        setattr(p_fit, name, float(v % (2*np.pi)))

    yhat = evaluate(t, p_fit, model, mix, compound, strict, harm)
    ss_res = float(np.sum((y - yhat)**2))
    ss_tot = float(np.sum((y - y.mean())**2))
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 0.0
    rmse = float(np.sqrt(ss_res/n))
    mae  = float(np.mean(np.abs(y - yhat)))

    # De-aliased band variance (interesting diagnostic): residual computed
    # against the pure lunar skeleton, after subtraction-filtering the data.
    if mix.semiannual or mix.evection:
        y_dem  = demixed_data(t, y, p_fit, model, mix, strict)
        yhat_dem = evaluate_demixed(t, p_fit, model, mix, compound, strict, harm)
        bv_dem = gf._band_variance(t, y_dem - yhat_dem)
    else:
        bv_dem = None

    diag = dict(
        r2=r2, rmse=rmse, mae=mae, ss_res=ss_res, ss_tot=ss_tot,
        n_samples=n,
        n_params=len(x0),
        param_names=list(names),
        nfev=int(res.nfev),
        cost=float(res.cost),
        status=int(res.status),
        residual_band_variance=gf._band_variance(t, y - yhat),
        residual_band_variance_demixed=bv_dem,
        mix_active=mix.label,
        compound_active=compound.label,
        strict_orbital=bool(strict.orbital),
        harmonics_active=harm.label,
    )

    if verbose:
        bits = []
        if mix.semiannual or mix.evection: bits.append(f"mix={mix.label}")
        if compound.lunar:                 bits.append(f"compound={compound.label}")
        if strict.orbital:                 bits.append("strict-orbital")
        if harm.level != "none":           bits.append(f"harmonics={harm.label}")
        suffix = (", " + ", ".join(bits)) if bits else ""
        print(f"[fit ]  closed-form ({model}{suffix}): R^2 = {r2:.6f}, "
              f"RMSE = {rmse:.4g}, "
              f"params = {len(x0)}, nfev = {res.nfev}")

    return p_fit, diag


# ============================================================================
# Output
# ============================================================================
def save_outputs(out_dir: str, t: np.ndarray, y: np.ndarray,
                 yhat: np.ndarray, p: ClosedParams, diag: dict,
                 boot: dict, model: str, time_unit: str,
                 mix: MixSpec | None = None,
                 compound: CompoundSpec | None = None,
                 strict: StrictSpec | None = None,
                 harm: HarmonicSpec | None = None):
    if mix is None:      mix = MixSpec()
    if compound is None: compound = CompoundSpec()
    if strict is None:   strict = StrictSpec()
    if harm is None:     harm = HarmonicSpec()
    os.makedirs(out_dir, exist_ok=True)

    bv = diag.get("residual_band_variance", {})
    var_total = float(np.var(y))
    var_resid = float(np.var(y - yhat))

    lines = [
        f"Closed-form lunisolar fit  ({model})",
        "=" * 60,
        f"  n samples           : {diag['n_samples']}",
        f"  closed-form params  : {diag['n_params']}",
        f"  R^2                 : {diag['r2']:.6f}",
        f"  RMSE                : {diag['rmse']:.6g}",
        f"  MAE                 : {diag['mae']:.6g}",
        "",
        "Comparison vs. linear-basis fit",
        "-" * 60,
        f"  linear-basis R^2    : {boot['linear_r2']:.6f}    "
        f"(unconstrained, ~80 amplitudes)",
        f"  closed-form  R^2    : {diag['r2']:.6f}    "
        f"({diag['n_params']} parameters)",
        f"  R^2 gap             : {boot['linear_r2'] - diag['r2']:+.6f}",
        "",
        "  The R^2 gap quantifies how much variance the strict closed",
        "  form gives up vs. the unconstrained lattice fit. Small gap =>",
        "  the closed form's locked amplitude relationships hold in the",
        "  data. Large gap => the data has structure (asymmetric side-",
        "  bands, even harmonics, independent forcings) that the closed",
        "  form cannot represent.",
        "",
        "Recovered closed-form parameters",
        "-" * 60,
        f"  a            = {p.a:.6g}        (carrier amplitude)",
        f"  c_0          = {p.c_0:+.6f}        (DC of inclination factor; "
        f"enables standalone anomalistic at 27.55 d)",
        f"  i_eff        = {p.i_eff:.6f} rad   "
        f"({np.rad2deg(p.i_eff):.3f} deg)   (effective inclination)",
        f"  e_eff        = {p.e_eff:.6f}        (effective eccentricity)",
        f"  m_N          = {p.m_N:.6f}        (nodal modulation depth)",
        f"  d_drac       = {p.d_drac:.4f} rad   "
        f"({np.rad2deg(p.d_drac):.2f} deg)",
        f"  d_ano        = {p.d_ano:.4f} rad   "
        f"({np.rad2deg(p.d_ano):.2f} deg)",
        f"  d_N          = {p.d_N:.4f} rad   "
        f"({np.rad2deg(p.d_N):.2f} deg)",
        f"  a_sid        = {p.a_sid:+.6g}       (sidereal block amp; "
        f"P_sid=27.3216 d. Modulated by (1+e_eff*cos th_a) so it also "
        f"spawns ano+sid at 13.719 d.)",
        f"  d_sid        = {p.d_sid:.4f} rad   "
        f"({np.rad2deg(p.d_sid):.2f} deg)",
        f"  b0           = {p.b0:.6g}",
    ]
    if model == "extended":
        lines += [
            "",
            "  Extended-model additions",
            f"  i2           = {p.i2:.6g}       (even-harmonic amplitude, "
            f"breaks Jacobi-Anger no-2*drac restriction)",
            f"  d_drac2      = {p.d_drac2:.4f} rad   "
            f"({np.rad2deg(p.d_drac2):.2f} deg)",
        ]
    if compound.lunar:
        lines += [
            "",
            "  Compound-tide additions  (additive; canonical periods locked)",
            f"  A_Mf         = {p.A_Mf:+.6g}        "
            f"@ T_Mf   = {T_Mf} d   (lunar fortnightly, K1-O1)",
            f"  phi_Mf       = {p.phi_Mf:.4f} rad   "
            f"({np.rad2deg(p.phi_Mf):.2f} deg)",
            f"  A_Mt         = {p.A_Mt:+.6g}        "
            f"@ T_Mt   = {T_Mt} d    (ter-mensal lunar)",
            f"  phi_Mt       = {p.phi_Mt:.4f} rad   "
            f"({np.rad2deg(p.phi_Mt):.2f} deg)",
            f"  A_MSqm       = {p.A_MSqm:+.6g}        "
            f"@ T_MSqm = {T_MSqm} d    (compound quarter-monthly)",
            f"  phi_MSqm     = {p.phi_MSqm:.4f} rad   "
            f"({np.rad2deg(p.phi_MSqm):.2f} deg)",
        ]
    if strict.orbital:
        lines += [
            "",
            "  --strict-orbital ACTIVE: i2 suppressed in extended model;",
            "  orbital module is the 11-param pure form (carrier with c_0",
            "  + sidereal block). Compound-tide module carries the",
            "  13.6 / 9 / 6.8 d compound-tide band variance.",
        ]
    if harm.level != "none":
        lines += [
            "",
            f"  Free-harmonic block (level={harm.label}): explicit additive ",
            f"  amplitudes on the lattice (breaking Bessel ratios).",
        ]
        for label, kd, kN, ka in _harm_lattice(harm.level):
            A = float(getattr(p, "A_" + label))
            phi = float(getattr(p, "phi_" + label))
            P = 2*np.pi/_harm_omega(kd, kN, ka) if _harm_omega(kd, kN, ka) != 0 else float('inf')
            lines += [
                f"  {label:<22s} A={A:+.6g}  phi={phi:.3f} rad  "
                f"P={P:.4f} d  (kd,kN,ka)=({kd:+d},{kN:+d},{ka:+d})"
            ]
    if mix.semiannual or mix.evection:
        lines += [
            "",
            "  Mixing-carrier additions  (multiplicative; subtraction-filtered)",
        ]
        if mix.semiannual:
            lines += [
                f"  alpha_sa     = {p.alpha_sa:+.6f}        "
                f"(semi-annual modulation depth)",
                f"  T_sa         = {p.T_sa:.4f} d         "
                f"(canonical 182.6211 d)",
                f"  phi_sa       = {p.phi_sa:.4f} rad   "
                f"({np.rad2deg(p.phi_sa):.2f} deg)",
            ]
        if mix.evection:
            lines += [
                f"  alpha_ev     = {p.alpha_ev:+.6f}        "
                f"(evection-aliased modulation depth)",
                f"  T_ev         = {p.T_ev:.4f} d           "
                f"(canonical 31.811 d)",
                f"  phi_ev       = {p.phi_ev:.4f} rad   "
                f"({np.rad2deg(p.phi_ev):.2f} deg)",
            ]
    lines += [
        "",
        "Residual band-variance check",
        "-" * 60,
        f"  total Var(y)              : {var_total:.4g}",
        f"  Var(residual)             : {var_resid:.4g}    "
        f"({100*var_resid/var_total:.2f}% of total)",
        f"    draconic     27.2122 d  : {bv.get('drac',0):.4g}",
        f"    sidereal     27.3216 d  : {bv.get('sidereal',0):.4g}",
        f"    sid_mirror   27.1036 d  : {bv.get('sid_mirror',0):.4g}",
        f"    anomalistic  27.5545 d  : {bv.get('anomalistic',0):.4g}",
        f"    nodal        6798.4 d   : {bv.get('nodal',0):.4g}",
        f"    after subtracting all 5 : {bv.get('after_all',0):.4g}",
    ]
    bv_dem = diag.get("residual_band_variance_demixed")
    if bv_dem is not None:
        lines += [
            "",
            "  De-aliased residual (data minus mixing envelope, vs lunar skeleton)",
            f"    draconic     27.2122 d  : {bv_dem.get('drac',0):.4g}",
            f"    sidereal     27.3216 d  : {bv_dem.get('sidereal',0):.4g}",
            f"    sid_mirror   27.1036 d  : {bv_dem.get('sid_mirror',0):.4g}",
            f"    anomalistic  27.5545 d  : {bv_dem.get('anomalistic',0):.4g}",
            f"    nodal        6798.4 d   : {bv_dem.get('nodal',0):.4g}",
            f"    after subtracting all 5 : {bv_dem.get('after_all',0):.4g}",
        ]
    lines += [
        "",
        "Locked periods (days)",
        f"  P_DRAC      = {gf.P_DRACONIC}",
        f"  P_ANO       = {gf.P_ANOMALISTIC}",
        f"  P_NODAL     = {gf.P_NODAL}",
        f"  P_SID  (=drac - N)  = {gf.P_SIDEREAL}        (emergent)",
        f"  P_SID' (=drac + N)  = {gf.P_SID_MIRROR:.4f}  (emergent)",
        "",
        f"  input time unit : {time_unit}",
    ]
    metrics_path = os.path.join(out_dir, "metrics.txt")
    with open(metrics_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[save]  {metrics_path}")
    print("\n".join(lines))

    # JSON
    json_path = os.path.join(out_dir, "params.json")
    with open(json_path, "w") as f:
        json.dump({
            "model": model,
            "mix":   asdict(mix),
            "compound": asdict(compound),
            "strict":   asdict(strict),
            "params": asdict(p),
            "diagnostics": diag,
            "bootstrap": boot,
            "periods": dict(P_DRACONIC=gf.P_DRACONIC,
                            P_ANOMALISTIC=gf.P_ANOMALISTIC,
                            P_NODAL=gf.P_NODAL,
                            P_SIDEREAL_EMERGENT=gf.P_SIDEREAL,
                            P_SID_MIRROR_EMERGENT=gf.P_SID_MIRROR),
        }, f, indent=2, default=lambda o: list(o))
    print(f"[save]  {json_path}")

    # Residuals
    res_path = os.path.join(out_dir, "residual.dat")
    np.savetxt(res_path,
               np.column_stack([t, y, yhat, y - yhat]),
               header="t_days  y  model  residual", comments="# ")
    print(f"[save]  {res_path}")

    # Plot
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
    ax.plot(t[:n_show], yhat[:n_show], lw=1.0, color="C2",
            label=f"closed-form model ({model})")
    ax.set_xlabel("t (days)")
    ax.set_ylabel("y")
    title = (f"Closed-form fit ({model})  -  R^2 = {diag['r2']:.4f}, "
             f"a = {p.a:.3g}, i_eff = {np.rad2deg(p.i_eff):.2f} deg, "
             f"m_N = {p.m_N:.3f}, e_eff = {p.e_eff:.3f}")
    title_bits = []
    if mix.semiannual:
        title_bits.append(f"a_sa={p.alpha_sa:+.3f}")
    if mix.evection:
        title_bits.append(f"a_ev={p.alpha_ev:+.3f}")
    if compound.lunar:
        title_bits.append(f"A_Mf={p.A_Mf:+.3g}")
    if title_bits:
        title += "  [" + ", ".join(title_bits) + "]"
    if strict.orbital:
        title += "  strict"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(t[:n_show], (y - yhat)[:n_show], lw=0.7, color="C0")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("t (days)")
    ax.set_ylabel("residual")
    ax.set_title("residual = data - model")
    ax.grid(alpha=0.3)

    ax = axes[2]
    resid = y - yhat
    if len(resid) > 16 and gf._HAVE_SCIPY:
        try:
            from scipy.signal import lombscargle
            t_f = t.astype(float)
            n = len(t_f)
            window = 0.5 - 0.5*np.cos(2*np.pi*np.arange(n)/max(n-1,1))
            x = (resid - resid.mean()) * window
            freqs = np.linspace(1.0/28.0, 1.0/26.0, 4000)
            angfreqs = 2*np.pi*freqs
            pgram = lombscargle(t_f, x.astype(float), angfreqs, normalize=True)
            periods = 1.0/freqs
            ax.plot(periods, pgram, lw=0.8)
            for P, name, c in [
                (gf.P_DRACONIC,    "drac",         "C3"),
                (gf.P_SIDEREAL,    "sid",          "C2"),
                (gf.P_SID_MIRROR,  "drac+N",       "C2"),
                (gf.P_ANOMALISTIC, "ano",          "C1"),
            ]:
                if periods.min() <= P <= periods.max():
                    ax.axvline(P, color=c, lw=0.6, ls="--", alpha=0.7)
                    ymax = ax.get_ylim()[1]
                    ax.text(P, ymax*0.95, name, rotation=90, va="top",
                            ha="right", fontsize=8, color=c)
            ax.set_xlabel("period (days)")
            ax.set_ylabel("normalized power")
            ax.set_title("residual spectrum")
            ax.grid(alpha=0.3)
        except Exception as exc:
            ax.text(0.5, 0.5, f"LS error: {exc}",
                    ha="center", va="center", transform=ax.transAxes)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, "fit.png")
    fig.savefig(fig_path, dpi=110)
    plt.close(fig)
    print(f"[save]  {fig_path}")


# ============================================================================
# CLI
# ============================================================================
def main():
    ap = argparse.ArgumentParser(
        description=("Constrained nonlinear closed-form lunisolar fitter, "
                     "warm-started from the linear-basis OLS solver."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("input", nargs="?", default=None,
                    help="whitespace-delimited 't y' or single-column y. "
                         "Omit for built-in synthetic test.")
    ap.add_argument("--time-unit", choices=["days", "years", "seconds"],
                    default="days")
    ap.add_argument("--out-dir", default="closed_out")
    ap.add_argument("--model", choices=["pure", "extended"],
                    default="extended",
                    help="'pure' = strict 8-param product form (mostly "
                         "diagnostic; cannot fit even harmonics). "
                         "'extended' (default) = adds explicit even-harmonic "
                         "and 3*sid pieces so the model can actually reach "
                         "the lattice points the data populates.")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--n-days", type=int, default=8000)
    ap.add_argument("--noise", type=float, default=0.005)
    ap.add_argument("--seed",  type=int, default=42)
    ap.add_argument("--mix-semiannual", action="store_true",
                    help="enable multiplicative semi-annual mixing carrier "
                         "(~182.6 d), with subtraction filter for residual")
    ap.add_argument("--mix-evection", action="store_true",
                    help="enable multiplicative evection-band mixing carrier "
                         "(~31.81 d), with subtraction filter for residual")
    ap.add_argument("--synthetic-mix", choices=["none", "semiannual",
                                                "evection", "both"],
                    default="none",
                    help="if --synthetic, optionally inject a known mixing "
                         "carrier so we can verify recovery")
    ap.add_argument("--compound-tides", choices=["none", "lunar"],
                    default="none",
                    help="enable additive compound-tide module: 'lunar' "
                         "adds Mf (13.661 d), Mt (9.133 d), MSqm (6.831 d) "
                         "with locked periods and free amplitude+phase "
                         "each (6 extra params)")
    ap.add_argument("--strict-orbital", action="store_true",
                    help="in extended model, suppress i2 so "
                         "the orbital module is the strict 8-param product. "
                         "Combine with --compound-tides lunar to test "
                         "whether the even-harmonic peaks are geometric "
                         "or compound-tide contamination.")
    ap.add_argument("--synthetic-compound", choices=["none", "lunar"],
                    default="none",
                    help="if --synthetic, optionally inject known compound-"
                         "tide signals (Mf, Mt, MSqm) so we can verify "
                         "recovery")
    ap.add_argument("--harmonics", choices=["none", "k3", "full"],
                    default="none",
                    help="add explicit free-amplitude harmonic block on the "
                         "(kd, kN, ka) lattice. 'k3' adds the 3*drac, 2*drac "
                         "and 2*ano amplitudes (7 amp/phase pairs = 14 "
                         "params), breaking the strict Bessel ratios. "
                         "'full' adds asymmetric carrier sidebands and the "
                         "largest cross-band gaps (12 amp/phase = 24 params)."
                    )
    args = ap.parse_args()

    mix      = MixSpec(semiannual=args.mix_semiannual,
                       evection=args.mix_evection)
    compound = CompoundSpec(lunar=(args.compound_tides == "lunar"))
    strict   = StrictSpec(orbital=args.strict_orbital)
    harm     = HarmonicSpec(level=args.harmonics)

    if args.synthetic or args.input is None:
        print(f"[data]  using synthetic signal (n_days={args.n_days}, "
              f"noise={args.noise})")
        t, y = gf.synthetic_signal(n_days=args.n_days,
                                   noise=args.noise, seed=args.seed)
        # Optionally inject mixing carriers into the synthetic signal
        # to verify the fitter can recover them.
        if args.synthetic_mix in ("semiannual", "both"):
            inj_sa = 1.0 + 0.18 * np.cos(2*np.pi*t/T_SEMIANNUAL_DEFAULT + 0.7)
            y = y * inj_sa
            print(f"[data]  injected semiannual mix: alpha=0.18, "
                  f"T={T_SEMIANNUAL_DEFAULT} d, phi=0.7 rad")
        if args.synthetic_mix in ("evection", "both"):
            inj_ev = 1.0 + 0.10 * np.cos(2*np.pi*t/T_EVECTION_DEFAULT + 1.3)
            y = y * inj_ev
            print(f"[data]  injected evection mix: alpha=0.10, "
                  f"T={T_EVECTION_DEFAULT} d, phi=1.3 rad")
        if args.synthetic_compound == "lunar":
            # Inject canonical compound tides to test recovery.
            cmp_inj = (
                0.0030 * np.cos(2*np.pi*t/T_Mf   + 0.4) +
                0.0015 * np.cos(2*np.pi*t/T_Mt   + 1.1) +
                0.0008 * np.cos(2*np.pi*t/T_MSqm + 2.2)
            )
            y = y + cmp_inj
            print(f"[data]  injected lunar compound tides: "
                  f"A_Mf=0.003, A_Mt=0.0015, A_MSqm=0.0008")
        time_unit = "days (synthetic)"
        os.makedirs(args.out_dir, exist_ok=True)
        np.savetxt(os.path.join(args.out_dir, "synthetic_input.dat"),
                   np.column_stack([t, y]),
                   header="t_days  y", comments="# ")
    else:
        t, y, time_unit = gf.read_input(args.input, args.time_unit)
        print(f"[data]  read {len(y)} samples from {args.input}  "
              f"(time_unit={time_unit})")
        if len(t) > 1:
            span = t.max() - t.min()
            print(f"        span = {span:.1f} d "
                  f"= {span/gf.DAYS_PER_YEAR:.2f} yr "
                  f"({span/gf.P_NODAL:.2f} nodal cycles)")

    p0, boot = bootstrap_from_linear(y, t, verbose=True)
    bits = []
    if mix.semiannual or mix.evection: bits.append(f"mix={mix.label}")
    if compound.lunar:                 bits.append(f"compound={compound.label}")
    if strict.orbital:                 bits.append("strict-orbital")
    if harm.level != "none":           bits.append(f"harmonics={harm.label}")
    suffix = (", " + ", ".join(bits)) if bits else ""
    print(f"[fit ]  refining closed-form ({args.model}{suffix}) with TRF...")
    p_fit, diag = fit_closed(y, t, p0, model=args.model,
                             mix=mix, compound=compound, strict=strict,
                             harm=harm, verbose=True)

    yhat = evaluate(t, p_fit, args.model, mix, compound, strict, harm)
    save_outputs(args.out_dir, t, y, yhat, p_fit, diag, boot,
                 args.model, time_unit,
                 mix=mix, compound=compound, strict=strict, harm=harm)


if __name__ == "__main__":
    main()
