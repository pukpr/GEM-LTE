import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 1. Load the data
df = pd.read_csv('NvsT.txt', sep='\s+')

# Extract columns
m_n = df['ModelN'].values
m_t = df['ModelT'].values
time_axis = np.arange(len(df))  # Common time axis (sample index)

# 2. Create the 3D Visualization
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the points as a scatter to see density
# Color is mapped to time to show the progression of the "orbit"
scatter = ax.scatter(m_n, m_t, time_axis, c=time_axis, cmap='viridis', s=10, alpha=0.6)

# Optionally plot the trajectory line to see the path between states
ax.plot(m_n, m_t, time_axis, color='gray', linewidth=0.5, alpha=0.3)

# 3. Labeling and Aesthetics
ax.set_xlabel('ModelN')
ax.set_ylabel('ModelT')
ax.set_zlabel('Time (Sample Index)')
ax.set_title('3D System Trajectory: ModelN vs ModelT over Time')

# Add a colorbar to indicate time progression
cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)
cbar.set_label('Time Progression')

# Adjust the viewing angle to emphasize the stratified layers
ax.view_init(elev=25, azim=45)

plt.tight_layout()
plt.show()

# 4. Secondary View: Top-Down "Orbit" Analysis
fig2 = plt.figure(figsize=(10, 8))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(m_n, m_t, time_axis, c=df['ForcingN'], cmap='plasma', s=8)
ax2.set_title('3D Trajectory Colored by Forcing Magnitude')
ax2.set_xlabel('ModelN')
ax2.set_ylabel('ModelT')
ax2.set_zlabel('Time')
ax2.view_init(elev=90, azim=-90) # Top-down perspective
plt.show()