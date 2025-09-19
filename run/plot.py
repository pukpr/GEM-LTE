import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.image as mpimg


identifier = sys.argv[1]
exquad = sys.argv[2]
start = sys.argv[3]
stop = sys.argv[4]
display = sys.argv[5]

# Read the PNG file (replace with your file name)
img = mpimg.imread('locs/' + identifier + '_loc.png')

# Load main data
df = pd.read_csv('lte_results.csv')
time = df.iloc[:, 0]
model = df.iloc[:, 1]
data = df.iloc[:, 2]
freq = df.iloc[:, 4]
model_psd = df.iloc[:, 5]
data_psd = df.iloc[:, 6]

# Read legend label
with open('lte_label.txt', 'r') as f:
    label_text = f.read()
with open('metrics.txt', 'r') as f:
    metrics = f.read()

# Read training points (expects two lines: X,Y)
training_data = np.loadtxt('training.txt', delimiter=',')
# Format: training_data[0] = [X1, Y1], training_data[1] = [X2, Y2]

fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=False, facecolor='#e6f2ff')

# Top chart: Time Series
axs[0].plot(time, model, linewidth=1, color='red', label='Model')
axs[0].plot(time, data, linewidth=1, color='blue', label='Data', alpha=0.7)
axs[0].set_title('Site#' + identifier + ' Time Series, biennial limit='+exquad)
axs[0].set_xlabel('Year')
axs[0].set_ylabel('Value')
axs[0].legend(loc='lower right', bbox_to_anchor=(0.98, 0.02))

# Draw thick dashed line for training points (top chart)
axs[0].plot(training_data[:, 0], training_data[:, 1], 'k--', linewidth=3, label='Training Segment')

# Add label in a box, lower left
props = dict(boxstyle='round', facecolor='white', alpha=0.8)
axs[0].text(0.02, 0.02, label_text, transform=axs[0].transAxes,
            fontsize=10, verticalalignment='bottom', horizontalalignment='left',
            bbox=props)

# Middle chart: Running Windowed Correlation
window = 50
if len(model) >= window:
    corrs = [np.corrcoef(model[i:i+window], data[i:i+window])[0, 1] for i in range(len(model)-window+1)]
    corr_time = time[window-1:]
    # Add an inset axes at [left, bottom, width, height] in axes coordinates (0-1)
    #inset_ax = axs[1].inset_axes([-0.1, 0.01, 0.5, 0.3])  # Adjust position and size as needed
    #inset_ax.imshow(img)
    #inset_ax.axis('off')  # Hide the axes for the inset
    axs[1].plot(corr_time, corrs, linewidth=3, color='green', zorder=1)
    axs[1].set_title(f'Running Windowed Correlation (window={window} months)')
    axs[1].set_xlabel('Year')
    axs[1].set_ylabel('Correlation Coefficient')
    # Draw thick dashed line for training points (middle chart)
    axs[1].plot(training_data[:, 0], training_data[:, 1], 'k--', linewidth=3, label='Training Segment', zorder=2)
    axs[1].imshow(img, extent=axs[1].get_xlim() + axs[1].get_ylim(), aspect='auto', alpha=0.25, zorder=0)
else:
    axs[1].text(0.5, 0.5, 'Not enough data for running correlation', ha='center', va='center')
    axs[1].set_axis_off()

# Bottom chart: Power Spectrum (log/log)
axs[2].loglog(freq, model_psd, linewidth=1, color='red', label='Model PSD')
axs[2].loglog(freq, data_psd, linewidth=1, color='blue', label='Data PSD', alpha=0.7)
axs[2].set_title('Power Spectrum')
axs[2].set_xlabel('Frequency (1/year)')
axs[2].set_ylabel('Power')
axs[2].legend(loc='upper right')
axs[2].text(0.02, 0.02, metrics, transform=axs[2].transAxes,
            fontsize=6, verticalalignment='bottom', horizontalalignment='left',
            bbox=props)


	    
# plt.rcParams['figure.facecolor'] = '#e6f2ff'

plt.tight_layout()

if display == '1':
    plt.show()
else:
    plt.savefig(identifier+'site'+start+'-'+stop+'.png', bbox_inches='tight')  # Save as PNG with tight bounding box

