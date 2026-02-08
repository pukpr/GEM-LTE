import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# -----------------------------
# Parameters (nondimensional)
# -----------------------------

omega0 = 1.0          # free nutation frequency (sets timescale)
zeta   = 1e-3         # weak damping

omega_m = 30.0        # "draconic" forcing frequency (off-resonant)
eps_m   = 2e-3        # forcing amplitude

Ts = 2*np.pi * 13.0   # sampling period (annual analogue)
eps_s = 3e-3          # impulse strength

tmax = 0.04 #4000.0         # total integration time
dt_sample = 0.5       # output sampling

# -----------------------------
# ODE system between impulses
# -----------------------------

def rhs(t, y):
    px, py, vx, vy = y

    fx = eps_m * np.cos(omega_m * t)
    fy = eps_m * np.sin(omega_m * t)

    dpx = vx
    dpy = vy
    dvx = -2*zeta*vx - omega0**2 * px + fx
    dvy = -2*zeta*vy - omega0**2 * py + fy

    return [dpx, dpy, dvx, dvy]

# -----------------------------
# Impulse event
# -----------------------------

def impulse_event(t, y):
    return np.sin(np.pi * t / Ts)

impulse_event.terminal = True
impulse_event.direction = 0

# -----------------------------
# Apply impulse map
# -----------------------------

def apply_impulse(y):
    px, py, vx, vy = y
    vx -= eps_s * omega0**2 * px
    vy -= eps_s * omega0**2 * py
    return np.array([px, py, vx, vy])

# -----------------------------
# Time integration loop
# -----------------------------

t = 0.0
y = np.array([1e-3, 0.0, 0.0, 0.0])  # small initial tilt

T_hist = []
P_hist = []

while t < tmax:
    sol = solve_ivp(
        rhs,
        (t, tmax),
        y,
        events=impulse_event,
        max_step=0.2,
        rtol=1e-9,
        atol=1e-9
    )

    # store solution
    T_hist.append(sol.t)
    P_hist.append(sol.y)

    # advance time
    t = sol.t[-1]
    y = sol.y[:, -1]

    # apply impulse if event triggered
    if sol.status == 1:
        y = apply_impulse(y)

# concatenate results
T = np.concatenate(T_hist)
P = np.hstack(P_hist)

px, py = P[0], P[1]

# -----------------------------
# Diagnostics
# -----------------------------

r = np.sqrt(px**2 + py**2)

# Remove initial transient
mask = T > 0.2 * tmax
T2 = T[mask]
r2 = r[mask]

# FFT
dt = np.mean(np.diff(T2))
freq = np.fft.rfftfreq(len(r2), dt)
spec = np.abs(np.fft.rfft(r2 - np.mean(r2)))

# -----------------------------
# Plots
# -----------------------------

plt.figure(figsize=(12,4))

plt.subplot(1,2,1)
plt.plot(T2, r2, lw=0.8)
plt.xlabel("Time")
plt.ylabel("|p|")
plt.title("Emergent Polar Motion Amplitude")

plt.subplot(1,2,2)
plt.semilogy(freq, spec)
plt.axvline(omega_m/(2*np.pi), color='r', ls='--', label="forcing")
plt.xlabel("Frequency")
plt.ylabel("Power")
plt.title("Spectrum")
plt.legend()

plt.tight_layout()
plt.show()
