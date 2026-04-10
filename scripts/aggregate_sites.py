import argparse
import urllib.request
from html.parser import HTMLParser

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from scipy.interpolate import UnivariateSpline

# -----------------------------
# Region filter helpers
# -----------------------------

class _PSMSLTableParser(HTMLParser):
    """Minimal HTML parser that extracts rows from the PSMSL station table."""
    def __init__(self):
        super().__init__()
        self.in_tbody = False
        self.in_tr = False
        self.in_td = False
        self.cells = []
        self.current_cell = ""
        self.rows = []

    def handle_starttag(self, tag, attrs):
        if tag == "tbody":
            self.in_tbody = True
        if self.in_tbody and tag == "tr":
            self.in_tr = True
            self.cells = []
        if self.in_tr and tag == "td":
            self.in_td = True
            self.current_cell = ""

    def handle_endtag(self, tag):
        if tag == "tbody":
            self.in_tbody = False
        if self.in_tbody and tag == "tr":
            if len(self.cells) >= 4:
                self.rows.append(self.cells[:])
            self.in_tr = False
        if self.in_tr and tag == "td":
            self.cells.append(self.current_cell.strip())
            self.in_td = False

    def handle_data(self, data):
        if self.in_td:
            self.current_cell += data


def _fetch_psmsl_stations():
    """Return list of (station_id, lat, lon) from the PSMSL RLR station table."""
    url = "https://psmsl.org/data/obtaining/index.php"
    with urllib.request.urlopen(url) as response:
        html = response.read().decode("utf-8")
    parser = _PSMSLTableParser()
    parser.feed(html)
    stations = []
    for row in parser.rows:
        try:
            stations.append((int(row[1]), float(row[2]), float(row[3])))
        except (ValueError, IndexError):
            pass
    return stations


_REGION_FILTERS = {
    "baltic-north-sea": None,
    "central-south-pacific": None,
    "north-atlantic": None,
    "australia": None,
    "west-coast-usa": None,
    "japan": None,
    "mediterranean": None,
    "southeast-asia": None,
    "hawaii": None,
    "uk-west-brittany": None,
    "east-coast-usa": None,
}


def get_region_station_ids(region):
    """Return a set of station IDs belonging to the named region."""
    if region not in _REGION_FILTERS:
        supported = ", ".join(_REGION_FILTERS)
        raise ValueError(f"Unknown region '{region}'. Supported: {supported}")

    print(f"Fetching PSMSL station list to filter for '{region}' region...")
    stations = _fetch_psmsl_stations()

    if region == "baltic-north-sea":
        def in_region(lat, lon):
            # North Sea: ~51–62°N, 5°W–9°E
            north_sea = (51 <= lat <= 62) and (-5 <= lon <= 9)
            # Baltic Sea: ~53–66°N, 9–30°E
            baltic_sea = (53 <= lat <= 66) and (9 <= lon <= 30)
            # Scandinavia (Norway, Sweden, Denmark, Finland): ~55–72°N, 4–32°E
            scandinavia = (55 <= lat <= 72) and (4 <= lon <= 32)
            return north_sea or baltic_sea or scandinavia

    elif region == "central-south-pacific":
        def in_region(lat, lon):
            # Pacific spans the dateline: western side (lon >= 110°E) and
            # eastern side (lon <= -65°W), latitude -60°S to 35°N.
            # This covers Australia, New Zealand, Pacific Island groups,
            # Hawaii, and the South/Central Pacific basin.
            in_lat = (-60 <= lat <= 35)
            in_lon = (lon >= 110) or (lon <= -65)
            return in_lat and in_lon

    elif region == "australia":
        def in_region(lat, lon):
            # Australian mainland + Tasmania: lat -44 to -10°S, lon 110–155°E.
            # Extended west to 96°E to include Cocos (Keeling) Islands (~97°E)
            # and east to 160°E to include Lord Howe Island (~159°E).
            return (-44 <= lat <= -10) and (96 <= lon <= 160)


    elif region == "north-atlantic":
        def in_region(lat, lon):
            # North Atlantic basin: equator to ~72°N, American coast to
            # western Europe/Africa (lon -100 to +5°E, stopping before
            # the Mediterranean).
            if not ((0 <= lat <= 72) and (-100 <= lon <= 5)):
                return False
            # Exclude St. Lawrence Seaway inland stations (Quebec City,
            # Trois-Rivières, Tadoussac, etc.) — seaway runs west of ~67°W
            # above 46°N.
            if lat > 46 and lon < -67:
                return False
            # Exclude Hudson Bay (not open Atlantic; e.g. Churchill, MB).
            if lat > 55 and lon < -75:
                return False
            return True

    elif region == "west-coast-usa":
        def in_region(lat, lon):
            # Continental US Pacific coast: San Diego (~32°N) to the
            # Canadian border (~49°N), longitude band -126 to -116°W.
            return (32 <= lat <= 55) and (-126 <= lon <= -116)

    elif region == "japan":
        def in_region(lat, lon):
            # Japanese archipelago including Ryukyu Islands (Okinawa/Ishigaki
            # ~24°N) through Hokkaido (~45.5°N).
            return (24 <= lat <= 46) and (122 <= lon <= 148)

    elif region == "mediterranean":
        def in_region(lat, lon):
            # Mediterranean + Black Sea basin (lon -5.5 to 42°E, lat 30-47.5°N).
            # Exclude Bay of Biscay / Atlantic Iberian coast: stations above
            # 42°N west of 0° are Atlantic, not Mediterranean.
            if not ((30 <= lat <= 47.5) and (-5.5 <= lon <= 42)):
                return False
            if lat > 42 and lon < 0:   # Bay of Biscay / Atlantic Spain/France
                return False
            return True

    elif region == "southeast-asia":
        def in_region(lat, lon):
            # Southeast Asia: Indochina, Malay Peninsula, Indonesian
            # archipelago, Philippines. Capped at 23°N to exclude southern
            # China coast and Japanese Ryukyu Islands (handled by 'japan').
            return (-10 <= lat <= 23) and (95 <= lon <= 135)

    elif region == "hawaii":
        def in_region(lat, lon):
            # Hawaiian Island chain: lat 18-23°N, lon -161 to -154°W.
            return (18 <= lat <= 23) and (-161 <= lon <= -154)

    elif region == "uk-west-brittany":
        def in_region(lat, lon):
            # Atlantic/Irish-Sea-facing coasts of the UK, Ireland, Isle of Man,
            # Channel Islands and Brittany (W France): lat 47-60°N, lon -11 to
            # -1.5°W.  Excludes the east-facing Scottish coast (Firth of Forth,
            # Moray Firth: Leith, Rosyth, Dunbar, Aberdeen, Buckie, Wick).
            if not ((47 <= lat <= 60) and (-11 <= lon <= -1.5)):
                return False
            if lat > 55 and lon > -4:   # east Scottish coast
                return False
            return True

    elif region == "east-coast-usa":
        def in_region(lat, lon):
            # Continental US Atlantic coast: Florida Keys (~24°N) to the
            # Maine/Canada border (~45°N), longitude band -82 to -65°W.
            return (24 <= lat <= 45) and (-82 <= lon <= -65)

    ids = {sid for sid, lat, lon in stations if in_region(lat, lon)}
    print(f"  Found {len(ids)} stations in '{region}' region.")
    return ids


# -----------------------------
# Regional average output
# -----------------------------

def write_regional_average(region_ids, region_name):
    """
    Load the full 1880-2026 dataset, apply the region filter, then process in
    six consecutive 25-year windows (1880-1905, 1905-1930, …, 2005-2030).

    Within each window every station is demeaned over that window to remove
    long-term trend and inter-station bias.  At each time step the average is
    sum(demeaned values) / number of contributing stations.  Time steps with
    no contributing stations are omitted from the output.

    The six sections are concatenated, dates are corrected from the PSMSL
    mid-month convention (year + (month-0.5)/12) to a start-of-month
    convention (year + (month-1)/12) so the series begins at exactly 1880.0,
    and finally a 12-point centred boxcar filter is applied to suppress the
    annual cycle.

    Output: <region_name>.dat, two space-separated columns (date, amplitude).
    """
    df_full = pd.read_csv(
        "filtered_psmsl_1880_2026.txt",
        delim_whitespace=True,
        header=None,
        names=["station", "date", "amp"],
    )
    df_full = df_full[df_full["amp"] != -99999]

    if region_ids is not None:
        df_full = df_full[df_full["station"].isin(region_ids)]
        if df_full.empty:
            print(f"No data for region '{region_name}' in filtered_psmsl_1880_2026.txt.")
            return

    # Six 25-year windows spanning 1880-2030
    windows = [(1880 + 25 * i, 1880 + 25 * (i + 1)) for i in range(6)]

    sections = []
    for t_start, t_end in windows:
        win = df_full[
            (df_full["date"] >= t_start) & (df_full["date"] < t_end)
        ].copy()

        if win.empty:
            continue

        # Demean each station within this window to remove its local bias
        win["amp_dm"] = (
            win["amp"] - win.groupby("station")["amp"].transform("mean")
        )

        # Average across all contributing stations at each time step
        avg = win.groupby("date")["amp_dm"].mean().reset_index()
        avg.columns = ["date", "avg"]
        sections.append(avg)

    if not sections:
        print("No data to average.")
        return

    result = pd.concat(sections).sort_values("date").reset_index(drop=True)

    # Correct dates: PSMSL mid-month  →  start-of-month
    # PSMSL: year + (month - 0.5)/12   →   year + (month - 1)/12
    # Equivalently: round(date * 12 - 0.5) / 12
    result["date"] = np.round(result["date"] * 12 - 0.5) / 12

    # 12-point boxcar (uniform) filter to suppress the annual cycle
    result["avg"] = result["avg"].rolling(window=12, center=True, min_periods=12).mean()
    result = result.dropna(subset=["avg"]).reset_index(drop=True)

    outfile = f"{region_name}.dat"
    with open(outfile, "w") as fh:
        for _, row in result.iterrows():
            fh.write(f"{row['date']:.4f}  {row['avg']:.4f}\n")

    print(f"Regional average written to '{outfile}' ({len(result)} time steps, "
          f"{result['date'].min():.4f} – {result['date'].max():.4f}).")


# -----------------------------
# CLI argument parsing
# -----------------------------
parser = argparse.ArgumentParser(
    description="Plot PSMSL sea-level time series with optional region filtering."
)
parser.add_argument(
    "--region",
    choices=["baltic-north-sea", "central-south-pacific", "north-atlantic",
             "australia", "west-coast-usa", "japan", "mediterranean",
             "southeast-asia", "hawaii", "uk-west-brittany", "east-coast-usa"],
    default=None,
    help="Only plot stations belonging to the specified region "
         "(fetches station coordinates from psmsl.org).",
)
parser.add_argument(
    "--start-date",
    type=float,
    default=2020.0,
    metavar="YEAR",
    help="Only include data on or after this decimal year (default: 2020.0).",
)
parser.add_argument(
    "--end-date",
    type=float,
    default=2030.0,
    metavar="YEAR",
    help="Only include data before this decimal year (default: 2030.0).",
)

parser.add_argument(
    "--align-signs",
    action="store_true",
    default=False,
    help="Flip the sign of each station's time series so its correlation with "
         "the ensemble is maximised. The first station in the data sets the "
         "reference polarity; subsequent stations are iteratively aligned to "
         "the running mean until convergence.",
)
parser.add_argument(
    "--average-output",
    action="store_true",
    default=False,
    help="Compute a regional average time series and write it to <region>.dat "
         "(or global.dat). Loads the full 1880-2026 dataset, processes it in "
         "six 25-year sections to remove long-term trend, then stitches and "
         "writes one row per contributing time step.",
)
args = parser.parse_args()

# -----------------------------
# Load data
# -----------------------------
# Expected columns:
# station_id   date_float   amplitude
df = pd.read_csv(
    "filtered_psmsl_1950_2026.txt",
    delim_whitespace=True,
    header=None,
    names=["station", "date", "amp"]
)

# Remove invalid values
df = df[df["amp"] != -99999]

# Apply region filter if requested
if args.region:
    region_ids = get_region_station_ids(args.region)
    df = df[df["station"].isin(region_ids)]
    if df.empty:
        print("No data found for the selected region. Exiting.")
        raise SystemExit(1)
else:
    region_ids = None

region_name = args.region if args.region else "global"

# Apply start-date / end-date filter
df = df[(df["date"] >= args.start_date) & (df["date"] < args.end_date)]
if df.empty:
    print(f"No data found between {args.start_date} and {args.end_date}. Exiting.")
    raise SystemExit(1)

# Write regional average if requested (uses independent 1880-2026 data load)
if args.average_output:
    write_regional_average(region_ids, region_name)

# -----------------------------
# Pre-compute processed series (demean + remove per-station annual cycle)
# -----------------------------
station_series = {}   # station_id -> (x_array, y_array)
for station, group in df.groupby("station"):
    group = group.sort_values("date")
    x = group["date"].values
    y = group["amp"].values - group["amp"].mean()

    # Fit annual sinusoid individually: amplitude and phase vary with location
    # (hemisphere flips sign; longitude shifts phase).
    X = np.column_stack([
        np.sin(2*np.pi*x),
        np.cos(2*np.pi*x)
    ])
    B, C = np.linalg.lstsq(X, y, rcond=None)[0]
    y = y - (B * np.sin(2*np.pi*x) + C * np.cos(2*np.pi*x))

    station_series[station] = (x, y)

# -----------------------------
# Optional sign alignment
# -----------------------------
def _align_signs(series_dict):
    """
    Iteratively flip each station's sign so it correlates positively with the
    current ensemble mean.  The first station fixes the reference polarity.
    Repeats until no further flips are needed (typically 1–2 passes).
    """
    stations = list(series_dict.keys())
    if len(stations) < 2:
        return series_dict

    # Round dates to 4 d.p. to allow index-based alignment across stations
    # (monthly decimal-year steps are ~0.083 apart so 4 d.p. is unambiguous).
    indexed = {
        sid: pd.Series(y, index=np.round(x, 4))
        for sid, (x, y) in series_dict.items()
    }

    signs = {sid: 1 for sid in stations}

    for iteration in range(10):  # safety cap; converges in ≤ 2 passes
        # Build ensemble mean from current signed series
        all_dates = sorted({d for sid in stations for d in indexed[sid].index})
        mean_series = pd.Series(0.0, index=all_dates)
        count_series = pd.Series(0, index=all_dates)
        for sid in stations:
            s = indexed[sid] * signs[sid]
            mean_series = mean_series.add(s, fill_value=0)
            count_series = count_series.add(pd.Series(1, index=s.index), fill_value=0)
        mean_series = mean_series / count_series.replace(0, np.nan)

        flips = 0
        for sid in stations:
            s = indexed[sid] * signs[sid]
            common = mean_series.dropna().index.intersection(s.index)
            if len(common) < 3:
                continue
            r = np.corrcoef(mean_series[common].values, s[common].values)[0, 1]
            if r < 0:
                signs[sid] *= -1
                flips += 1

        if flips == 0:
            print(f"  Sign alignment converged after {iteration + 1} iteration(s).")
            break
    else:
        print("  Sign alignment: reached iteration cap.")

    total_flipped = sum(1 for sid in stations if signs[sid] == -1)
    print(f"  {total_flipped} of {len(stations)} station(s) had their sign flipped.")

    return {
        sid: (series_dict[sid][0], series_dict[sid][1] * signs[sid])
        for sid in stations
    }, {sid for sid in stations if signs[sid] == -1}


if args.align_signs:
    print("Aligning signs across stations...")
    station_series, flipped_ids = _align_signs(station_series)
else:
    flipped_ids = set()

# -----------------------------
# Plot setup
# -----------------------------
plt.figure(figsize=(12, 8))

# Loop through each station
for station, (x, y) in station_series.items():

    # Plot faint lines (red if sign was flipped, steelblue otherwise)
    color = "red" if station in flipped_ids else "steelblue"
    plt.plot(x, y, color=color, alpha=0.25, linewidth=1)

    # Identify outlier stations
    if np.max(np.abs(y)) > 400:

        # Index of the absolute max anomaly in the raw demeaned data
        idx = np.argmax(np.abs(y))

        # x-value of the max anomaly
        x_peak = x[idx]

        # Evaluate spline at that location
        y_peak = y[idx]

        k = station % 6
        # Vertical alignment cycling
        va_options = ["bottom", "top", "center", "baseline", "center_baseline",
                  "bottom", "top", "center", "baseline", "center_baseline"]
        va = va_options[k]

        # Compute variation index
        k = station % 10

        # Horizontal alignment cycling
        ha_options = ["right", "left", "center", "left", "right",
                     "center", "left", "right", "center", "left"]
        ha = ha_options[k]

        k = 10*(random.random()-0.5)
        k1 = 10*(random.random()-0.5)

        # Offset label position (varies slightly with k)
        x_offset = 0.1 + 0.1 * (k - 2)
        y_offset = 50 + 30 * (k1 - 2)
 
        # Label slightly offset from the peak
        plt.annotate(
            f"{station}",
            xy=(x_peak, y_peak),
            # xytext=(x_peak + 0.1, y_peak + 50),   # offset label
            xytext=(x_peak + x_offset, y_peak + y_offset),
            arrowprops=dict(
                arrowstyle="-",
                color="black",
                lw=0.1
            ),
            fontsize=6,
            ha=ha,
            va=va
        )


# -----------------------------
# Final formatting
# -----------------------------
n_stations = len(station_series)
if args.region:
    region_label = f"{args.region}  ({n_stations} stations)"
else:
    region_label = f"global  ({n_stations} stations)"
plt.title(f"Sea-Level Anomaly Time Series — {region_label}\nDemeaned, annual cycle removed per station")
plt.xlabel("Year")
plt.ylabel("Amplitude (mm)")
plt.grid(True, alpha=0.2)

plt.show()


