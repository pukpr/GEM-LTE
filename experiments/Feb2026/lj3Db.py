import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_3d_time_evolution(file_path):
    # 1. Load Data
    try:
        df = pd.read_csv(file_path, sep='\s+')
    except FileNotFoundError:
        print("File not found.")
        return

    # 2. Construct Common Forcing Axis (Depth)
    # We still use the average forcing for the Z-axis to unfold the loops
    df['Forcing_Common'] = (df['ForcingN'] + df['ForcingT']) / 2.0

    # 3. 3D Visualization
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Data mapping
    xs = df['ModelN']
    ys = df['ModelT']
    zs = df['Forcing_Common']
    
    # COLOR MAPPING: Use the DataFrame index (Row Number)
    # This shows the progression from the first data point to the last.
    colors = df.index 

    # Plot scatter
    # 'viridis' or 'plasma' are good for time sequences (dark to light)
    scatter = ax.scatter(xs, ys, zs, c=colors, cmap='turbo', s=6, alpha=0.7, depthshade=True)

    # Labeling
    ax.set_xlabel('ModelN Response')
    ax.set_ylabel('ModelT Response')
    ax.set_zlabel('Common Forcing')
    ax.set_title('3D Lissajous Trajectory\nColor = Time (Sequence Order)')

    # Add a color bar to indicate the time progression
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.7)
    cbar.set_label('Sequence Index (Start -> End)')

    # Set initial view angle
    # Elev=30 looks down slightly to see the loops
    # Azim=-60 rotates to see the spread between models
    ax.view_init(elev=30, azim=-60)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_3d_time_evolution('NvsT.txt')