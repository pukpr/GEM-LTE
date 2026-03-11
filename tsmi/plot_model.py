import csv
import json
import matplotlib.pyplot as plt
import numpy as np

t = []
I = []
I_model = []

# Read params for title
try:
    with open('ts.json', 'r') as f:
        params = json.load(f)
        zeta = params['damping']['zeta']
        omega0 = params['damping']['omega0']
        annual_phase = params.get('annual_phase', 0.0)
        ic_A = params.get('initial_condition', {}).get('A', 0.0)
        ic_B = params.get('initial_condition', {}).get('B', 0.0)
except Exception as e:
    print(f"Warning: Could not read ts.json: {e}")
    zeta = 0
    omega0 = 0
    annual_phase = 0

try:
    with open('model_fit.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t.append(float(row['t']))
                I.append(float(row['I']))
                I_model.append(float(row['I_model']))
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

plt.figure(figsize=(12, 6))
plt.plot(t, I, label='Observed (I)', color='blue', alpha=0.6, linewidth=1)
plt.plot(t, I_model, label='Model (I_model)', color='red', alpha=0.6, linewidth=1)

plt.xlabel('Time')
plt.ylabel('Value')
plt.title(f'Time Series: Observed vs Model (r = {r_val:.4f})\n$\\zeta$={zeta:.3f}, $\\omega_0$={omega0:.3f}, $\\phi$={annual_phase:.3f}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

output_file = 'model_fit.png'
plt.savefig(output_file, dpi=150)
print(f'Plot saved to {output_file}')
