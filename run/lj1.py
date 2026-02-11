import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# 1. Load Data
# Note: The file contains irregular spacing/tabs, sep='\s+' handles this.
df = pd.read_csv('NvsT.txt', sep='\s+')

# Extract values
fn, mn = df['ForcingN'].values, df['ModelN'].values
ft, mt = df['ForcingT'].values, df['ModelT'].values

# 2. Linear Stretching & Compressing Optimization
def objective(params):
    scale, offset = params
    # Apply stretching/compression and offset to the T-axis forcing
    ft_adj = scale * ft + offset
    
    # Define common range for comparison
    common_min = max(fn.min(), ft_adj.min())
    common_max = min(fn.max(), ft_adj.max())
    
    if common_min >= common_max: return 1e6 # Invalid range penalty
    
    common_x = np.linspace(common_min, common_max, len(df))
    
    # Interpolate both onto the adjusted common forcing axis
    m_n_interp = interp1d(fn, mn, bounds_error=False, fill_value="extrapolate")(common_x)
    m_t_interp = interp1d(ft_adj, mt, bounds_error=False, fill_value="extrapolate")(common_x)
    
    # Calculate correlation (we minimize negative correlation to maximize)
    correlation = np.corrcoef(m_n_interp, m_t_interp)[0, 1]
    return -abs(correlation)

# Initial guess: scale=1, offset=0
initial_guess = [1.0, 0.0]
res = minimize(objective, initial_guess, method='Nelder-Mead')
best_scale, best_offset = res.x

print(f"Optimal Scale: {best_scale:.4f}")
print(f"Optimal Offset: {best_offset:.4f}")
print(f"Maximized Correlation: {-res.fun:.4f}")

# Generate Optimized Data
ft_optimized = best_scale * ft + best_offset
common_x = np.linspace(max(fn.min(), ft_optimized.min()), min(fn.max(), ft_optimized.max()), len(df))
mn_opt = interp1d(fn, mn)(common_x)
mt_opt = interp1d(ft_optimized, mt)(common_x)

# 3. Dynamic Time Warping (DTW) - Simplified Implementation
# This finds a non-linear path to align the sequences
from scipy.spatial.distance import cdist

def simple_dtw(s1, s2):
    n, m = len(s1), len(s2)
    dtw_matrix = np.zeros((n+1, m+1))
    dtw_matrix[1:, 0] = np.inf
    dtw_matrix[0, 1:] = np.inf
    cost_matrix = cdist(s1.reshape(-1,1), s2.reshape(-1,1), metric='euclidean')
    for i in range(1, n+1):
        for j in range(1, m+1):
            cost = cost_matrix[i-1, j-1]
            dtw_matrix[i, j] = cost + min(dtw_matrix[i-1, j], dtw_matrix[i, j-1], dtw_matrix[i-1, j-1])
    return dtw_matrix

# Note: For very large datasets, use the 'fastdtw' library instead.
# For this script, we use the optimized linear results for visualization.

# 4. Visualization
plt.figure(figsize=(15, 5))

# Plot 1: Original vs Optimized Lissajous
plt.subplot(1, 3, 1)
plt.plot(mn, mt, 'r--', alpha=0.3, label='Original')
plt.plot(mn_opt, mt_opt, 'b-', label='Optimized (Linear)')
plt.title('Lissajous: Original vs Optimized')
plt.xlabel('ModelN')
plt.ylabel('ModelT')
plt.legend()

# Plot 2: Forcing Comparison
plt.subplot(1, 3, 2)
plt.plot(fn, label='ForcingN', alpha=0.7)
plt.plot(ft_optimized, label='ForcingT (Stretched)', alpha=0.7)
plt.title('Aligned Forcing Axes')
plt.legend()

# Plot 3: Aligned Responses
plt.subplot(1, 3, 3)
plt.plot(common_x, mn_opt, label='ModelN')
plt.plot(common_x, mt_opt, label='ModelT')
plt.title('Time-Synchronized Responses')
plt.legend()

plt.tight_layout()
plt.show()