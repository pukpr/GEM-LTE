import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import matplotlib.image as mpimg
import warnings

identifier = sys.argv[1]
marker = sys.argv[2]
start = sys.argv[3]
stop = sys.argv[4]
display = sys.argv[5]

# Get the directory where the script itself is located
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Build the path relative to the script
img_path = os.path.join(BASE_DIR, 'locs', f"{identifier}_loc.png")

# Read the PNG file (replace with your file name)
img = mpimg.imread(img_path)
# img = mpimg.imread('locs/' + identifier + '_loc.png')

# Load main data
df = pd.read_csv('lte_results.csv')
time = df.iloc[:, 0]
model = df.iloc[:, 1]
data = df.iloc[:, 2]
forcing = df.iloc[:, 3]
freq = df.iloc[:, 4]
model_psd = df.iloc[:, 5]
data_psd = df.iloc[:, 6]

mean_forcing_path = os.path.join(BASE_DIR, 'mean_forcing.dat')
mean_forcing = None
if os.path.isfile(mean_forcing_path):
    try:
        mean_data = np.loadtxt(mean_forcing_path, ndmin=2)
        if mean_data.shape[1] >= 2:
            mean_forcing = mean_data[:, :2]
    except Exception:
        mean_forcing = None

label_text = identifier

#with open('metrics.txt', 'r') as f:
#    metrics = f.read()

start_time = float(start) # training_data[0,0]
stop_time = float(stop) # training_data[1,0]

fig, axs = plt.subplots(3, 2, figsize=(16, 12), gridspec_kw={'width_ratios': [2.5, 1]}, sharex=False, facecolor='#e6f2ff')

# Top chart: Time Series
axs[0,0].plot(time, model, linewidth=1, color='red', label='Model')
axs[0,0].plot(time, data, linewidth=1, color='blue', label='Data', alpha=0.7)
axs[0,0].set_title('Site#' + identifier + ' Time Series, cross-validation:'+marker)
axs[0,0].set_xlabel('Year')
axs[0,0].set_ylabel('Value')
axs[0,0].legend(loc='lower right', bbox_to_anchor=(0.98, 0.02))

# Draw thick dashed line for training points (top chart)
#axs[0,0].plot(training_data[:, 0], training_data[:, 1], 'k--', linewidth=3, label='Training Segment')
axs[0,0].plot([start_time, stop_time], [0.0, 0.0], 'k--', linewidth=3, label='Training Segment')



# Middle chart: Forcing
axs[1,0].plot(time, forcing, linewidth=1, color='red')
if mean_forcing is not None:
    mean_mask = (mean_forcing[:, 0] >= time.min()) & (mean_forcing[:, 0] <= time.max())
    if np.any(mean_mask):
        axs[1,0].plot(
            mean_forcing[mean_mask, 0],
            mean_forcing[mean_mask, 1],
            'k--',
            linewidth=0.6,
            label='Mean forcing',
        )
        axs[1,0].legend(loc='upper left')
axs[1,0].set_title('Latent/hidden forcing layer')
axs[1,0].set_xlabel('Year')
axs[1,0].set_ylabel('Value')


# Bottom chart: Running Windowed Correlation
window = 50
if len(model) >= window:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        corrs = [np.corrcoef(model[i:i+window], data[i:i+window])[0, 1] for i in range(len(model)-window+1)]
    corr_time = time[window-1:]
    # Add an inset axes at [left, bottom, width, height] in axes coordinates (0-1)
    #inset_ax = axs[1].inset_axes([-0.1, 0.01, 0.5, 0.3])  # Adjust position and size as needed
    #inset_ax.imshow(img)
    #inset_ax.axis('off')  # Hide the axes for the inset
    axs[2,0].plot(corr_time, corrs, linewidth=3, color='green', zorder=1)
    axs[2,0].set_title(f'Running Windowed Correlation (window={window} months)')
    axs[2,0].set_xlabel('Year')
    axs[2,0].axvline(start_time, ls='--', label="forcing")
    axs[2,0].axvline(stop_time, ls='--', label="forcing")
    axs[2,0].set_ylabel('Correlation Coefficient')
    # Draw thick dashed line for training points (middle chart)
    axs[2,0].plot([start_time, stop_time], [0.0, 0.0], 'k--', linewidth=3, label='Training Segment', zorder=2)
    axs[2,0].imshow(img, extent=axs[2,0].get_xlim() + axs[2,0].get_ylim(), aspect='auto', alpha=0.25, zorder=0)
else:
    axs[2,0].text(0.5, 0.5, 'Not enough data for running correlation', ha='center', va='center')
    axs[2,0].set_axis_off()

# Add label in a box, lower left
props = dict(boxstyle='round', facecolor='yellow', alpha=0.8)
axs[2,0].text(0.02, 0.02, label_text, transform=axs[2,0].transAxes,
            fontsize=24, verticalalignment='bottom', horizontalalignment='left',
            bbox=props)


# Middle chart: Regression
cc = np.corrcoef(model, data)[0,1]

# Define your interval
# Create a boolean mask
mask = (time >= start_time) & (time <= stop_time)

# Subset both time and Y
time_subset = time[mask]
model_subset = model[mask]  # Slicing columns (time steps)
data_subset = data[mask]  # Slicing columns (time steps)
cc_valid = np.corrcoef(model_subset, data_subset)[0,1]

not_mask = (time < start_time) | (time > stop_time)
time_train = time[not_mask]
model_train = model[not_mask]  # Slicing columns (time steps)
data_train = data[not_mask]  # Slicing columns (time steps)
cc_train = np.corrcoef(model_train, data_train)[0,1]


axs[0,1].plot(model, data, linestyle="None", marker='.', color='green', label='training=' + f"{cc_train:.4f}" )
axs[0,1].plot(model_subset, data_subset, linestyle="None", marker='.', color='red', label='validation=' + f"{cc_valid:.4f}")
axs[0,1].set_title('Regression CC=' + f"{cc:.4f}")
axs[0,1].set_ylabel('Model')
axs[0,1].set_xlabel('Data')
axs[0,1].legend(loc='upper left')


# Middle chart: LTE Forcing Modulation
axs[1,1].plot(model, forcing, linestyle="None", marker='_', color='red', label='Model')
axs[1,1].plot(data, forcing,  linestyle="None", marker='_', color='blue', label='Data')
axs[1,1].set_title('LTE Modulation')
axs[1,1].set_ylabel('Latent forcing level')
axs[1,1].set_xlabel('Level Modulation')
axs[1,1].legend(loc='upper left')


# Bottom right chart: Power Spectrum (log/log)
axs[2,1].loglog(freq, model_psd, linewidth=1, color='red', label='Model PSD')
axs[2,1].loglog(freq, data_psd, linewidth=1, color='blue', label='Data PSD', alpha=0.7)
axs[2,1].set_title('Power Spectrum, modulation on latent')
axs[2,1].set_xlabel('Frequency (per level)')
axs[2,1].set_ylabel('Power')
axs[2,1].set_xlim(left=1.0)
axs[2,1].legend(loc='lower left')

#axs[2,1].text(0.02, 0.02, metrics, transform=axs[2,1].transAxes,
#            fontsize=6, verticalalignment='bottom', horizontalalignment='left',
#            bbox=props)



	    
# plt.rcParams['figure.facecolor'] = '#e6f2ff'

plt.tight_layout()

if display == '1':
    plt.show()
else:
    plt.savefig(identifier+'site'+start+'-'+stop+'.png', bbox_inches='tight')  # Save as PNG with tight bounding box
