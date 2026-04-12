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
    "south-america-pacific": None,
    "north-atlantic": None,
    "australia": None,
    "west-coast-usa": None,
    "japan": None,
    "mediterranean": None,
    "southeast-asia": None,
    "hawaii": None,
    "uk-west-brittany": None,
    "east-coast-usa": None,
    "gulf-of-mexico": None,
    "new-zealand": None,
    "canada-maritimes": None,
    "alaska-se": None,
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

    elif region == "south-america-pacific":
        def in_region(lat, lon):
            # Pacific span of the South America coast
            # eastern side (lon <= -65°W), latitude -60°S to 35°N.
            in_lat = (-60 <= lat <= 15)
            in_lon = (-100 <= lon <= -65)
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

    elif region == "gulf-of-mexico":
        def in_region(lat, lon):
            # Gulf of Mexico: Yucatan (~18°N) to the northern coast (~31°N),
            # Texas/Mexico coast (-98°W) to the Florida peninsula (-80°W).
            # The -80°W eastern boundary excludes Atlantic-facing Florida stations.
            return (18 <= lat <= 31) and (-98 <= lon <= -80)

    elif region == "new-zealand":
        def in_region(lat, lon):
            # North Island to southern South Island + Stewart Island: ~47–34°S,
            # 166–178°E.  Excludes Chatham Islands (~176°W) and Australian coast.
            return (-47 <= lat <= -34) and (166 <= lon <= 178)

    elif region == "canada-maritimes":
        def in_region(lat, lon):
            # Nova Scotia, New Brunswick (Bay of Fundy + Gulf of St. Lawrence),
            # Prince Edward Island, and coastal Newfoundland/Labrador.
            # Lat 43–52°N, lon -67 to -52°W keeps coastal Gulf/Atlantic stations
            # and excludes St. Lawrence River inland sites (lon < -67°W).
            return (43 <= lat <= 52) and (-67 <= lon <= -52)

    elif region == "alaska-se":
        def in_region(lat, lon):
            # Southeast Alaska panhandle (Ketchikan to Yakutat) plus
            # Prince Rupert / Dixon Entrance area of northern BC.
            # Lat 54–61°N, lon -141 to -129°W.
            return (54 <= lat <= 61) and (-141 <= lon <= -129)

    ids = {sid for sid, lat, lon in stations if in_region(lat, lon)}
    print(f"  Found {len(ids)} stations in '{region}' region.")
    return ids


def get_region_station_locs(region):
    """Return a dict {station_id: (lat, lon)} for all stations in the region."""
    stations = _fetch_psmsl_stations()
    if region is None:
        return {sid: (lat, lon) for sid, lat, lon in stations}
    region_ids = get_region_station_ids(region)
    return {sid: (lat, lon) for sid, lat, lon in stations if sid in region_ids}


def save_region_map(region_name, region_ids, active_ids):
    """
    Plot all stations in the region as grey dots; those that contributed data
    within the selected date interval are drawn as filled red dots.
    Saves to <region_name>_map.png.

    Parameters
    ----------
    region_name : str
        Used for the plot title and output filename.
    region_ids : set or None
        All station IDs in the region (None → global).
    active_ids : set
        Station IDs that have data in the filtered time window.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    all_stations = _fetch_psmsl_stations()
    if region_ids is not None:
        region_stations = [(sid, lat, lon) for sid, lat, lon in all_stations
                           if sid in region_ids]
    else:
        region_stations = all_stations

    lats_all = [lat for _, lat, _ in region_stations]
    lons_all = [lon for _, _, lon in region_stations]
    lats_act = [lat for sid, lat, lon in region_stations if sid in active_ids]
    lons_act = [lon for sid, lat, lon in region_stations if sid in active_ids]

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines(linewidth=0.6)
    ax.add_feature(cfeature.LAND, facecolor="lightgrey", edgecolor="black",
                   linewidth=0.4)
    ax.add_feature(cfeature.OCEAN, facecolor="white")

    # All region stations: grey background dots
    if lats_all:
        ax.plot(lons_all, lats_all, "o", color="brown", markersize=3,
                transform=ccrs.PlateCarree(), label="in region (no data in window)",
                zorder=3)

    # Stations active in the date window: red dots on top
    if lats_act:
        ax.plot(lons_act, lats_act, "o", color="red", markersize=4,
                transform=ccrs.PlateCarree(), label="active in date window",
                zorder=4)

    title = region_name if region_name != "global" else "Global"
    plt.title(f"{title}  —  {len(lats_act)} active / {len(lats_all)} total stations")
    ax.legend(loc="lower left", fontsize=8)
    plt.tight_layout()

    outfile = f"{region_name}_map.png"
    plt.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Station map saved to '{outfile}'.")


# -----------------------------
# Regional average output
# -----------------------------

def write_regional_average(region_ids, region_name, detrend=False):
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

        if detrend:
            def _detrend_station(grp):
                t = grp["date"].values
                if len(t) < 2:
                    return grp["amp_dm"]
                coeffs = np.polyfit(t, grp["amp_dm"].values, 1)
                return pd.Series(
                    grp["amp_dm"].values - np.polyval(coeffs, t),
                    index=grp.index,
                )
            win["amp_dm"] = win.groupby("station", group_keys=False).apply(_detrend_station)

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


def write_regional_average_all(region_ids, region_name, detrend=False):
    """
    Like write_regional_average but processes the entire time span in one pass
    (no 25-year windowing).  Each station is demeaned over its full record,
    then all stations are averaged at each time step, dates are corrected to
    start-of-month, and a 12-point centred boxcar filter is applied.

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

    # Demean each station over its entire record
    df_full = df_full.copy()
    df_full["amp_dm"] = (
        df_full["amp"] - df_full.groupby("station")["amp"].transform("mean")
    )

    if detrend:
        def _detrend_station(grp):
            t = grp["date"].values
            if len(t) < 2:
                return grp["amp_dm"]
            coeffs = np.polyfit(t, grp["amp_dm"].values, 1)
            return pd.Series(
                grp["amp_dm"].values - np.polyval(coeffs, t),
                index=grp.index,
            )
        df_full["amp_dm"] = df_full.groupby("station", group_keys=False).apply(_detrend_station)

    # Average across all contributing stations at each time step
    result = df_full.groupby("date")["amp_dm"].mean().reset_index()
    result.columns = ["date", "avg"]
    result = result.sort_values("date").reset_index(drop=True)

    # Correct dates: PSMSL mid-month  →  start-of-month
    result["date"] = np.round(result["date"] * 12 - 0.5) / 12

    # 12-point centred boxcar filter to suppress the annual cycle
    result["avg"] = result["avg"].rolling(window=12, center=True, min_periods=12).mean()
    result = result.dropna(subset=["avg"]).reset_index(drop=True)

    outfile = f"{region_name}.dat"
    with open(outfile, "w") as fh:
        for _, row in result.iterrows():
            fh.write(f"{row['date']:.4f}  {row['avg']:.4f}\n")

    print(f"Regional average (full span) written to '{outfile}' ({len(result)} time steps, "
          f"{result['date'].min():.4f} – {result['date'].max():.4f}).")


# -----------------------------
# CLI argument parsing
# -----------------------------
parser = argparse.ArgumentParser(
    description="Plot PSMSL sea-level time series with optional region filtering."
)
parser.add_argument(
    "--region",
    choices=["baltic-north-sea", "south-america-pacific", "north-atlantic",
             "australia", "west-coast-usa", "japan", "mediterranean",
             "southeast-asia", "hawaii", "uk-west-brittany", "east-coast-usa",
             "gulf-of-mexico", "new-zealand", "canada-maritimes", "alaska-se"],
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
parser.add_argument(
    "--detrend",
    action="store_true",
    default=False,
    help="After all filtering, remove a least-squares linear trend from the "
         "averaged output (applies to --average-output and --average-all).",
)
parser.add_argument(
    "--average-all",
    action="store_true",
    default=False,
    help="Like --average-output but demeans each station over its full record "
         "(no 25-year windowing). Writes <region>.dat after applying a "
         "12-point centred boxcar filter.",
)
parser.add_argument(
    "--save-map",
    action="store_true",
    default=False,
    help="Save a PNG map of all stations in the region. Grey dots = all stations "
         "in the region; red dots = stations with data in the selected date window. "
         "Output: <region>_map.png (or global_map.png).",
)

args = parser.parse_args()

# -----------------------------
# Load data
# -----------------------------
# Expected columns:
# station_id   date_float   amplitude
df = pd.read_csv(
    "filtered_psmsl_1880_2026.txt",
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

# Save station map if requested (active_ids = stations present after all filters)
if args.save_map:
    active_ids = set(df["station"].unique())
    save_region_map(region_name, region_ids, active_ids)

# Write regional average if requested (uses independent 1880-2026 data load)
if args.average_output:
    write_regional_average(region_ids, region_name, detrend=args.detrend)
if args.average_all:
    write_regional_average_all(region_ids, region_name, detrend=args.detrend)

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


