#!/usr/bin/env python3
"""comb_sines_10000yr.py

Creates a comb delta impulse (1/year) convolved with sum of five equal-amplitude sine waves,
simulates for 10,000 years at 1 sample/day, computes FFT of the result, and plots amplitude-squared.

Saves:
 - /mnt/data/comb_sines_power_10000yr.npz  (freqs, power)
 - /mnt/data/comb_sines_spectrum.png      (full-range plot)
 - /mnt/data/comb_sines_spectrum_zoom.png (zoom ~27-day region)

Usage: python3 /mnt/data/comb_sines_10000yr.py
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time
import os

# Parameters
YEARS = 10000
DAYS_PER_YEAR = 365.25
SAMPLES_PER_DAY = 1          # samples per day
dt_days = 1.0 / SAMPLES_PER_DAY

# Sinewave periods in days (as requested)
periods_days = np.array([27.2722, 27.2122/2.0, 27.3216, 27.3216/2.0, 27.5545])
amplitude_each = 1.0

# Derived sizes
total_days = int(round(YEARS * DAYS_PER_YEAR))
N = total_days * SAMPLES_PER_DAY
print(f"Total samples: N = {N:,} (approx {N/1e6:.3f} million)")

# Time axis (in days)
t = np.arange(N) * dt_days

# Build summed sine signal
s = np.zeros(N, dtype=np.float64)
for P in periods_days:
    s += amplitude_each * np.sin(2 * np.pi * t / P)

# Build comb of delta impulses at integer-year positions
comb = np.zeros(N, dtype=np.float64)
year_indices = np.round(np.arange(YEARS) * DAYS_PER_YEAR).astype(int)
year_indices = year_indices[year_indices < N]
comb[year_indices] = 1.0

# Convolve via multiplication in frequency domain (efficient)
t0 = time()
S_freq = np.fft.rfft(s)
Comb_freq = np.fft.rfft(comb)
R_freq = S_freq * Comb_freq
r = np.fft.irfft(R_freq, n=N)
t1 = time()
print(f"Convolution done in {t1-t0:.2f} s")

# FFT of result and amplitude-squared (power)
t0 = time()
R_freq2 = np.fft.rfft(r)
power = np.abs(R_freq2) ** 2
freqs_cpd = np.fft.rfftfreq(N, d=dt_days)   # cycles per day
freqs_cpy = freqs_cpd * DAYS_PER_YEAR       # cycles per year
t1 = time()
print(f"FFT computed in {t1-t0:.2f} s")

# Save spectrum
out_dir = "/tmp/fft"
os.makedirs(out_dir, exist_ok=True)
np.savez_compressed(os.path.join(out_dir, "comb_sines_power_10000yr.npz"),
                    freqs_cpy=freqs_cpy, power=power)
print("Saved power spectrum to", os.path.join(out_dir, "comb_sines_power_10000yr.npz"))

# Plot full-range spectrum (cycles per year)
plt.figure(figsize=(10,6))
plt.plot(freqs_cpy, power, linewidth=0.8)
plt.xlabel('Frequency (cycles per year)')
plt.ylabel('Amplitude squared (power)')
plt.title(f'Amplitude-squared spectrum — comb (1/yr) convolved with 5 sines — {YEARS:,} years')
plt.xlim(0, 30)
plt.ylim(bottom=0)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
full_png = os.path.join(out_dir, "comb_sines_spectrum.png")
plt.savefig(full_png, dpi=200)
plt.close()
print("Saved full-spectrum plot to", full_png)

# Plot zoom around ~27-day harmonics
peak_center = DAYS_PER_YEAR / 27.3  # approx cycles/year for ~27.3 days
lo = max(0, peak_center - 5)
hi = peak_center + 5
plt.figure(figsize=(10,5))
plt.plot(freqs_cpy, power, linewidth=0.8)
plt.xlabel('Frequency (cycles per year)')
plt.ylabel('Amplitude squared (power)')
plt.title(f'Zoom around ~27-day harmonics (approx {peak_center:.2f} cpy)')
plt.xlim(lo, hi)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
zoom_png = os.path.join(out_dir, "comb_sines_spectrum_zoom.png")
plt.savefig(zoom_png, dpi=200)
plt.close()
print("Saved zoom-spectrum plot to", zoom_png)

print('Done.')
