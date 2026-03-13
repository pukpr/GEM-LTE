import csv
import json
import matplotlib.pyplot as plt
import numpy as np

# Data lists
t = []
I = []
I_model = []
Manifold = []
Manifold_Empirical = []

try:
    with open('model_fit.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t.append(float(row['t']))
                I.append(float(row['I']))
                I_model.append(float(row['I_model']))
                if 'Manifold' in row:
                    Manifold.append(float(row['Manifold']))
                if 'Manifold_Empirical' in row:
                    try:
                        val = float(row['Manifold_Empirical'])
                        Manifold_Empirical.append(val)
                    except:
                        pass
            except ValueError:
                continue
except FileNotFoundError:
    print("Error: model_fit.csv not found. Run the Julia script first.")
    exit(1)

if not t:
    print("Error: No data found in model_fit.csv")
    exit(1)

# Calculate Pearson correlation coefficient
corr_matrix = np.corrcoef(I, I_model)
r_val = corr_matrix[0, 1]
print(f"Correlation coefficient (r): {r_val:.4f}")

# Create figure with 3 subplots (stacked)
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12), height_ratios=[1, 1, 1])

# Plot 1: Modulation Scatter (Manifold vs Output)
if Manifold:
    ax1.scatter(Manifold, I, label='Observed (I)', color='blue', s=10, alpha=0.5, marker='.')
    ax1.scatter(Manifold, I_model, label='Model (I_model)', color='red', s=10, alpha=0.5, marker='.')
    ax1.set_xlabel('Latent Manifold Value (X)')
    ax1.set_ylabel('Output Value (Y)')
    ax1.set_title('Modulation Function: Latent Manifold vs Output')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
else:
    ax1.text(0.5, 0.5, 'Manifold data not available', ha='center', va='center')

# Plot 2: Time Series
ax2.plot(t, I, label='Observed (I)', color='blue', alpha=0.6, linewidth=1)
ax2.plot(t, I_model, label='Model (I_model)', color='red', alpha=0.6, linewidth=1)
ax2.set_xlabel('Time')
ax2.set_ylabel('Value')
ax2.set_title(f'Time Series: Observed vs Model (r = {r_val:.4f})')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Time vs Latent Manifold
if Manifold:
    ax3.plot(t, Manifold, label='Parametric Manifold (Optimized)', color='green', linewidth=1)
    
    if Manifold_Empirical and len(Manifold_Empirical) == len(t):
        ax3.plot(t, Manifold_Empirical, label='Empirical Manifold (11force.dat)', color='gray', linestyle='--', linewidth=1, alpha=0.5)
        try:
             m_corr = np.corrcoef(Manifold, Manifold_Empirical)[0, 1]
             ax3.text(0.02, 0.95, f'Manifold Correlation: {m_corr:.4f}', transform=ax3.transAxes, fontsize=10, verticalalignment='top')
        except:
             pass
    
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Manifold Value')
    ax3.set_title('Latent Manifold Dynamics over Time')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
else:
    ax3.text(0.5, 0.5, 'Manifold data not available', ha='center', va='center')

plt.tight_layout()

output_file = 'model_fit.png'
plt.savefig(output_file, dpi=150)
print(f'Plot saved to {output_file}')
