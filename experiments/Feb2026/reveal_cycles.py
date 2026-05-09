import numpy as np
import pandas as pd
import glob
import os
import json
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import detrend

def cosine_model(t, A, f, p, c):
    """Synchronized sinusoidal model with phase relative to start of window."""
    return A * np.cos(2 * np.pi * f * t + p) + c

def analyze_directory(filepath):
    """Processes a results file and its accompanying parameter file."""
    dirpath = os.path.dirname(filepath) or "."
    p_file = os.path.join(dirpath, "lt.exe.p")
    
    try:
        # 1. Load Data (Col 2: Y, Col 4: X)
        df = pd.read_csv(filepath, usecols=[1, 3], header=None, engine='c')
        df.columns = ['Y', 'X']
        df = df.sort_values(by='X')
        cleaned = df.groupby('X')['Y'].mean().reset_index()
        x, y = cleaned['X'].values, cleaned['Y'].values
        span = x.max() - x.min()

        # 2. Spectral Initial Guess (FFT)
        # Resampled to 4096 points to handle non-uniformity and improve FFT resolution
        x_uni = np.linspace(x.min(), x.max(), 4096)
        y_uni = np.interp(x_uni, x, detrend(y))
        yf = np.fft.rfft(y_uni)
        xf = np.fft.rfftfreq(len(x_uni), x_uni[1] - x_uni[0])
        
        # Focus on the 3-10 cycle window for the "Principal Component"
        mask = (xf >= 3/span) & (xf <= 10/span)
        f_guess = xf[mask][np.argmax(np.abs(yf[mask]))] if np.any(mask) else xf[np.argmax(np.abs(yf[1:])) + 1]

        # 3. Sinusoidal Regression (Phase-Locked Sync)
        p0 = [np.std(y)*1.414, f_guess, 0, np.mean(y)]
        popt, _ = curve_fit(cosine_model, x, y, p0=p0)
        fit_A, fit_f, fit_p, fit_c = popt

        # Polarity Correction: ensure A is positive and move inversion info to phase
        if fit_A < 0: fit_A, fit_p = -fit_A, fit_p + np.pi
        fit_p = (fit_p + np.pi) % (2 * np.pi) - np.pi
        
        f_fit = abs(fit_f)
        period = 1 / f_fit
        cycles = span * f_fit

        # 4. Uncertainty Delta Calculation
        # Compare analyzed frequency against 'ltep' entries in the JSON parameter file
        u_delta = 0.0
        f_closest = f_fit
        if os.path.exists(p_file):
            with open(p_file, 'r') as f_json:
                p_data = json.load(f_json)
            
            # The 'ltep' field provides the set of expected/possible frequencies
            # ltep_vals = p_data.get('ltep', [])
            ltep_abs = [abs(v) for v in p_data.get('ltep', [])]
            ltep_vals = ltep_abs + [ltep_abs[-1] * h for h in p_data.get('harm', [])] if ltep_abs else []            
            if ltep_vals:
                f_closest = min(ltep_vals, key=lambda val: abs(val - f_fit))
                ref_period = 1.0/f_closest
                u_delta = abs(1.0/f_fit - 1.0/f_closest)

        # 5. Zero-Crossing Calculation
        t_hr = np.linspace(x.min(), x.max(), 5000)
        y_fit_hr = cosine_model(t_hr, fit_A, fit_f, fit_p, 0)
        zc_count = len(np.where(np.diff(np.sign(y_fit_hr)))[0])

        # 6. Output Result
        print(f"{dirpath}, {period:.6f}, {fit_p:.4f}, {cycles:.3f}, {zc_count}, {u_delta:.6f}")

        # 7. Visualization with Uncertainty Indicator
        plt.figure(figsize=(10, 5))
        plt.plot(x, y, 'k.', alpha=0.15, markersize=1, label='Raw Data')
        
        # Primary fit curve
        y_plot = cosine_model(t_hr, fit_A, fit_f, fit_p, fit_c)
        plt.plot(t_hr, y_plot, 'r-', linewidth=2, label=f'Sync Fit ({cycles:.2f} cyc)')
        
        # Uncertainty band/reference if applicable
        if u_delta > 0:
            y_ref = cosine_model(t_hr, fit_A, f_closest, fit_p, fit_c)
            plt.plot(t_hr, y_ref, 'b--', alpha=0.5, label=f'Ref λ ({ref_period:.4f})')
            plt.fill_between(t_hr, y_plot, y_ref, color='gray', alpha=0.15, label='Freq Uncertainty')

        # plt.title(f"Principal Sync: {dirpath}\nCycles: {cycles:.2f} | ZC: {zc_count} | Delta: {u_delta:.4f}")
        plt.title(f"Index: {dirpath}\nCycles: {cycles:.2f} | λ: {period:.4f} | Φ: {fit_p:0.4f}")
        plt.xlabel("X (Column 4)")
        plt.ylabel("Y (Column 2)")
        plt.legend(loc='upper right', fontsize='small')
        plt.grid(True, alpha=0.2)
        
        # Filename based on directory path
        plot_name = dirpath.replace(os.sep, '_').strip('.') + "_principal_fit.png"
        plt.savefig(plot_name or "root_fit.png")
        plt.close()

    except Exception as e:
        print(f"{dirpath}, ERROR: {str(e).splitlines()[0]}, 0, 0, 0, 0")

if __name__ == "__main__":
    # Required glob level: "*/lte_results.csv"
    files = sorted(glob.glob("*/lte_results.csv"))
    
    # Fallback to current directory if glob is empty
    if not files and os.path.exists("lte_results.csv"):
        files = ["lte_results.csv"]

    if not files:
        sys.exit(1)
        
    print("DirPath, Period, Phase_Rad, Cycles_Fractional, Zero_Crossings, Uncertainty_Delta")
    for f in files:
        analyze_directory(f)

