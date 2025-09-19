import sys
import numpy as np

def read_data(file):
    times = []
    values = []
    for line in file:
        if not line.strip():
            continue
        parts = line.strip().split(';')
        if len(parts) < 2:
            continue
        t = float(parts[0])
        v = float(parts[1])
        if v == -99999:
            continue
        times.append(t)
        values.append(v)
    return np.array(times), np.array(values)

def boxcar_filter(values, window_size=13):
    """Apply a boxcar (moving average) filter, handling endpoints with smaller windows."""
    n = len(values)
    filtered = np.zeros(n)
    half = window_size // 2
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        filtered[i] = np.mean(values[start:end])
    return filtered

def remove_linear_trend(times, values):
    """Remove linear trend using linear regression."""
    order = int(sys.argv[2])
    p = np.polyfit(times, values, order)  # 1=linear, 2=quadratic
    trend = np.polyval(p, times)
    detrended = values - trend
    return detrended

def main():
    if len(sys.argv) < 3:
        print("Usage: python boxcar_detrend_timeseries.py <data_file.txt>", file=sys.stderr)
        sys.exit(1)
    filename = sys.argv[1]
    with open(filename, 'r') as f:
        times, values = read_data(f)

    # Read from stdin
    # times, values = read_data(sys.stdin)
    if len(times) < 13:
        print("Not enough data points for 13-point boxcar filter.", file=sys.stderr)
        sys.exit(1)
    # Apply boxcar filter
    filtered = boxcar_filter(values, window_size=13)
    # Remove linear trend
    result = remove_linear_trend(times, filtered)
    # Output: time value (space-delimited)
    for t, v in zip(times, result):
        print(f"{t:.6f} {v:.6f}")

if __name__ == "__main__":
    main()
