import numpy as np
import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import detrend

def principal_sine_model(t, A, f, p, c):
    """The Principal Component model: a synchronized sine wave."""
    # We use (t - t_start) so the phase 'p' is relative to the start of the data.
    return A * np.sin(2 * np.pi * f * (t - t[0]) + p) + c

def analyze_file(filepath):
    try:
        # 1. Load Data (Column 2 -> Y, Column 4 -> X)
        df = pd.read_csv(filepath, usecols=[1, 3], header=None, engine='c')
        df.columns = ['Y', 'X']
        
        # 2. Sort by Column 4 and handle duplicate X entries
        df = df.sort_values(by='X')
        cleaned = df.groupby('X')['Y'].mean().reset_index()
        x, y = cleaned['X'].values, cleaned['Y'].values
        span = x.max() - x.min()

        # 3. Principal Frequency Detection (FFT)
        # Resample to uniform grid to ensure FFT accuracy
        x_uni = np.linspace(x.min(), x.max(), 4096)
        y_uni = np.interp(x_uni, x, detrend(y))
        
        yf = np.fft.rfft(y_uni)
        xf = np.fft.rfftfreq(len(x_uni), x_uni[1] - x_uni[0])
        
        # Target the 3-10 cycle frequency range
        mask = (xf >= 3/span) & (xf <= 10/span)
        if not np.any(mask):
            # Fallback if range is not present
            f_guess = xf[np.argmax(np.abs(yf[1:])) + 1]
        else:
            f_guess = xf[mask][np.argmax(np.abs(yf[mask]))]

        # 4. Sinusoidal Regression (The "Sync" Phase)
        # Refines the FFT guess into a perfect fit of the primary oscillation
        try:
            p0 = [np.std(y)*1.414, f_guess, 0, np.mean(y)]
            popt, _ = curve_fit(principal_sine_model, x, y, p0=p0)
            fit_A, fit_f, fit_p, fit_c = popt
        except:
            # Fallback to FFT guess if regression fails to converge
            fit_A, fit_f, fit_p, fit_c = np.std(y)*1.414, f_guess, 0, np.mean(y)

        # 5. Extract Parameters
        period = 1 / abs(fit_f)
        precise_cycles = span / period
        
        # Calculate Zero Crossings of the Principal Component (The cross-check)
        # We generate a high-res version of the fit to count crossings accurately
        t_highres = np.linspace(x.min(), x.max(), 5000)
        y_fit_centered = principal_sine_model(t_highres, fit_A, fit_f, fit_p, 0)
        zero_crossings = len(np.where(np.diff(np.sign(y_fit_centered)))[0])

        # 6. Output Result
        # Format: File, Period, Phase, Cycles, Zero_Crossings
        print(f"{filepath}, {period:.6f}, {fit_p:.4f}, {precise_cycles:.3f}, {zero_crossings}")

        # 7. Verification Plot
        plt.figure(figsize=(10, 5))
        plt.plot(x, y, 'k.', alpha=0.15, markersize=1, label='Raw Data (Sorted)')
        
        # Plot the Principal Component
        y_final_fit = principal_sine_model(t_highres, fit_A, fit_f, fit_p, fit_c)
        plt.plot(t_highres, y_final_fit, 'r-', linewidth=2, label='Principal Component')
        
        # Zero-crossing reference line
        plt.axhline(fit_c, color='blue', linestyle='--', alpha=0.5, label='Mean/Zero-X Line')
        
        plt.title(f"Principal Component: {filepath}\nCycles: {precise_cycles:.2f} | Crossings: {zero_crossings}")
        plt.xlabel("Column 4 (Sorted X)")
        plt.ylabel("Column 2 (Signal Y)")
        plt.legend(loc='upper right')
        plt.grid(True, alpha=0.2)
        
        # Save PNG
        plot_path = filepath.replace('.csv', '_principal.png').replace('/', '_')
        plt.savefig(plot_path)
        plt.close()

    except Exception as e:
        print(f"{filepath}, ERROR: {str(e).splitlines()[0]}, 0, 0, 0")

if __name__ == "__main__":
    # Traversal logic
    files = sorted(glob.glob("*/lte_results.csv", recursive=True))
    
    if not files:
        print("No files found.")
        sys.exit(1)
        
    print("File, Period, Phase_Rad, Cycles_Fractional, Zero_Crossings")
    for f in files:
        analyze_file(f)

