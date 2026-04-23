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
    return A * np.sin(2 * np.pi * f * (t - t[0]) + p) + c

def analyze_file(filepath):
    try:
        # Get the directory path only
        dirpath = os.path.dirname(filepath)
        
        # 1. Load Data (Column 2 -> Y, Column 4 -> X)
        df = pd.read_csv(filepath, usecols=[1, 3], header=None, engine='c')
        df.columns = ['Y', 'X']
        
        # 2. Sort by Column 4 and handle duplicate X entries
        df = df.sort_values(by='X')
        cleaned = df.groupby('X')['Y'].mean().reset_index()
        x, y = cleaned['X'].values, cleaned['Y'].values
        span = x.max() - x.min()

        # 3. Principal Frequency Detection (FFT)
        x_uni = np.linspace(x.min(), x.max(), 4096)
        y_uni = np.interp(x_uni, x, detrend(y))
        yf = np.fft.rfft(y_uni)
        xf = np.fft.rfftfreq(len(x_uni), x_uni[1] - x_uni[0])
        
        mask = (xf >= 3/span) & (xf <= 10/span)
        f_guess = xf[mask][np.argmax(np.abs(yf[mask]))] if np.any(mask) else xf[np.argmax(np.abs(yf[1:])) + 1]

        # 4. Sinusoidal Regression
        p0 = [np.std(y)*1.414, f_guess, 0, np.mean(y)]
        popt, _ = curve_fit(principal_sine_model, x, y, p0=p0)
        fit_A, fit_f, fit_p, fit_c = popt

        # 5. Extract Parameters
        period = 1 / abs(fit_f)
        precise_cycles = span / period
        
        t_highres = np.linspace(x.min(), x.max(), 5000)
        y_fit_centered = principal_sine_model(t_highres, fit_A, fit_f, fit_p, 0)
        zero_crossings = len(np.where(np.diff(np.sign(y_fit_centered)))[0])

        # 6. Output Result using dirpath
        print(f"{dirpath}, {period:.6f}, {fit_p:.4f}, {precise_cycles:.3f}, {zero_crossings}")

        # 7. Verification Plot
        plt.figure(figsize=(10, 5))
        plt.plot(x, y, 'k.', alpha=0.15, markersize=1, label='Raw Data')
        y_final_fit = principal_sine_model(t_highres, fit_A, fit_f, fit_p, fit_c)
        plt.plot(t_highres, y_final_fit, 'r-', linewidth=2, label='Principal Component')
        plt.axhline(fit_c, color='blue', linestyle='--', alpha=0.5)
        plt.title(f"Index: {dirpath}\nCycles: {precise_cycles:.2f} | λ: {period:.4f} | Φ: {fit_p:0.4f}")
        plt.grid(True, alpha=0.2)
        
        # Save PNG in the local directory with a name based on the path
        plot_name = dirpath + '/' + dirpath.replace(os.sep, '_') + "_fit.png" if dirpath else "root_fit.png"
        plt.savefig(plot_name)
        plt.close()

    except Exception as e:
        print(f"{os.path.dirname(filepath)}, ERROR: {str(e).splitlines()[0]}, 0, 0, 0")

if __name__ == "__main__":
    files = sorted(glob.glob("*/lte_results.csv", recursive=True))
    if not files:
        sys.exit(1)
        
    print("DirPath, Period, Phase_Rad, Cycles_Fractional, Zero_Crossings")
    for f in files:
        analyze_file(f)
        
