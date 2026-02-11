import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# 1. Load the data
# The file uses tabs as delimiters; 'sep=\s+' handles multiple spaces/tabs
file_path = 'NvsT.txt'
df = pd.read_csv(file_path, sep='\s+')

# 2. Extract columns
# Note: Ensure column names match the file header exactly
f_n = df['ForcingN'].values
m_n = df['ModelN'].values
f_t = df['ForcingT'].values
m_t = df['ModelT'].values

# 3. Define a common forcing range for interpolation
# We use the intersection of the ranges to avoid extrapolation errors
f_min = max(f_n.min(), f_t.min())
f_max = min(f_n.max(), f_t.max())
common_forcing = np.linspace(f_min, f_max, len(df))

# 4. Perform Interpolation
# This synchronizes ModelN and ModelT to the same forcing values
interp_n = interp1d(f_n, m_n, kind='linear', bounds_error=False, fill_value="extrapolate")
interp_t = interp1d(f_t, m_t, kind='linear', bounds_error=False, fill_value="extrapolate")

model_n_synced = interp_n(common_forcing)
model_t_synced = interp_t(common_forcing)

# 5. Visualization
plt.figure(figsize=(12, 5))

# Plot the synchronized time series (Forcing vs Models)
plt.subplot(1, 2, 1)
plt.plot(common_forcing, model_n_synced, label='ModelN (Synced)', alpha=0.7)
plt.plot(common_forcing, model_t_synced, label='ModelT (Synced)', alpha=0.7)
plt.xlabel('Common Forcing')
plt.ylabel('Model Response')
plt.title('Synchronized Responses')
plt.legend()
plt.grid(True)

# Plot the Lissajous Pattern (ModelN vs ModelT)
plt.subplot(1, 2, 2)
plt.plot(model_n_synced, model_t_synced, color='purple', lw=1)
plt.xlabel('ModelN')
plt.ylabel('ModelT')
plt.title('Lissajous Pattern (ModelN vs ModelT)')
plt.grid(True)

plt.tight_layout()
plt.show()

# Optional: Output the synchronized data to a new file
# synced_df = pd.DataFrame({
#     'Forcing_Common': common_forcing,
#     'ModelN_Synced': model_n_synced,
#     'ModelT_Synced': model_t_synced
# })
# synced_df.to_csv('Synchronized_Lissajous_Data.csv', index=False)
