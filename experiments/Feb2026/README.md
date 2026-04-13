# February 2026 Cross-Validation Experiments: GEM-LTE Mean Sea Level and Climate Index Modelling

## Executive Summary

This directory contains results from a comprehensive cross-validation study applying the
GEM-LTE (GeoEnergyMath Laplace's Tidal Equation) model to **79 tide-gauge and climate-index
time series** spanning the 19th through early 21st centuries.  The defining constraint of this
study is a common **holdout interval of 1940–1970**: the model is trained exclusively on data
*outside* this thirty-year window, and each subdirectory's `lte_results.csv` and
`*site1940-1970.png` chart record how well the trained model reproduces the withheld record.

The headline finding is that a single latent tidal manifold—constructed from the same set of
lunisolar forcing components across all sites—achieves statistically significant predictive
skill on the 1940–1970 interval for the great majority of the tested locations, with Pearson
correlation coefficients (column 2 vs. column 3 of `lte_results.csv`) ranging from r ≈ 0.72
at the best-performing Baltic tide gauges to r ≈ 0.12 at the most challenging Atlantic
stations.  Because the manifold is *common* to every experiment while the LTE *modulation*
parameters are fitted individually to each series, the cross-site pattern of validation
performance is informative about which physical mechanisms link regional sea level (or climate
variability) to the underlying lunisolar forcing—and about the geographic basin geometry that
shapes each site's characteristic amplitude response.

---

## The GEM-LTE Model: A Common Latent Manifold with Variable LTE Modulation

The model that generates every result in this directory is built in two stages, implemented
in `lt.exe` (the compiled Ada solver):

### 1 · The Latent Manifold

A scalar signal **M(t)** is constructed by numerically integrating a tidal forcing composed of
sinusoids at the principal lunisolar components (listed in the LPAP table of each `lt.exe.par`
file: 27.32 d, 27.21 d, 13.66 d, 13.63 d, 9.08 d, 6.82 d, ~1095 d, ~1616 d, ~6793 d, etc.).
The forcing is multiplied by an annual Gaussian comb—a train of narrow pulses, one per year,
whose width (`delB`) and offset (`delA ≈ 4.56`) are fixed across all tide-gauge experiments—
before being integrated over time.  The integration accumulates phase, producing a
slowly-varying continuous manifold that encodes the long-period envelope of the lunisolar
beating pattern. The superposition of the underlying long-period lunar tridal factors looks like the following -- note that this considers all (apart from QBO30, which is a special wavenumber=0 geometry) sites as varying colors:

<img width="871" height="308" alt="image" src="https://github.com/user-attachments/assets/f40ec7ac-ec0f-4f21-b7d6-ebc198ce505f" />

The bar chart of the strengths of the tidal factors is shown as the following with the spread indicated as varying colors:

<img width="499" height="536" alt="image" src="https://github.com/user-attachments/assets/7f975360-9db7-45c2-8046-a70c92b76784" />

The beat between the 27.212d day and 27.3216 day factor lead to an 18.6 year envolope, while the beat bewteen the strongest 27.2122 day and the 27.5545 day anomalistic factor leads to a weaker 6-year modulation. All the other lunar tidal factors reinforce these periods (with another 8.85 year modulation that is hard to distinguish from the 9.3 year harmonic).


A critical architectural feature is that the **background acceleration `bg = 1.4145 × 10⁻⁴`
yr⁻¹** is identical across all 77 experiments that share the same `lt.exe.par` template; only
the two QBO-band and two NINO-east experiments deviate slightly.  This universality means that
every site sees the same low-frequency secular ramp embedded in the manifold.

### 2 · LTE Modulation (Laplace's Tidal Equation)

The manifold M(t) is passed through a nonlinear LTE transformation—the nonlinear Laplace's
Tidal Equation solution whose mathematical basis is laid out in *Mathematical Geoenergy*
(Wiley/AGU, 2019; TeX source at https://github.com/pukpr/gem).  This transformation is
parameterised by a small set of `ltep` scale parameters and `harm` harmonics, optimised
independently for each site.

**The critical conceptual point is that the LTE modulation operates on the instantaneous
*amplitude level* of M(t), not on the time axis.**  The sinusoidal terms in the LTE
transformation take the form

    f(M) = Σᵢ sin( M / ltepᵢ + φᵢ )

where the argument of each sinusoid is the manifold *value* M divided by the scale parameter
ltepᵢ—not a frequency multiplied by time.  As a consequence, `ltep` values must **not** be
read as temporal periods or cycle durations.  They are **resonant wavenumber analogues** of
the fluidic bounding container: they describe the amplitude scale at which the enclosed-basin
geometry creates standing-wave modes in the fluid, in the same spirit as the spatial
wavenumbers that appear in Laplace's Tidal Equation solutions for a bounded fluid domain.
The numeric values of `ltep` happen to be expressed in year-equivalent units (matching the
units of M), but their physical interpretation is geometric—they encode how the depth
structure, coastline, and bathymetry of a particular ocean basin or enclosed sea project the
manifold amplitude onto observable sea-level (or climate-index) variance.

Examining the `lt.exe.par` files across all subdirectories reveals two recurring
site-specific wavenumber parameters:

- **ltep₁ ≈ ±23** (ranging roughly 22–24; polarity varies) — The dominant low-wavenumber
  resonant mode of the Baltic–North Sea enclosed basin system.  The vast majority of tide
  gauges around the Baltic and North Seas converge on a value near **23.2**.  The polarity
  (sign) reflects whether the LTE waveform is in-phase or anti-phase with the manifold at
  the particular basin geometry of each site, not a temporal phase difference.

- **ltep₂ ≈ 9.4–9.5** — A secondary resonant mode appearing as the dominant second
  modulation for virtually every Baltic and North Sea station, consistent with a higher
  spatial-harmonic of the enclosed-basin standing-wave response.

Beyond these two, a **third ltep = 18.68** is shared nearly universally as a fixed
secondary slot, accompanied by a set of long-scale slots (100, 200, 400, 600, 70, 150)
representing oceanic and climate variability at sub-Milankovitch amplitude scales, and a
short slot (0.5) for the semiannual.

Climate indices depart from this template more dramatically, reflecting their distinct
geographic basin settings: NINO4 carries ltep₁ = 12.1 (western equatorial Pacific warm-pool
geometry); the AMO uses ltep₁ = −7.2 and ltep₂ = −7.2 (a strongly symmetric two-mode
pattern in the North Atlantic basin); and NINO3 / NINO1+2 carry ltep₁ ≈ 94, indicating that
the eastern-Pacific coastal boundary responds at a very different spatial amplitude scale than
the western/central Pacific.

The interplay between the universal manifold and the site-specific LTE transformation is the
core of the model's generative power: the same underlying tide-driven forcing drives different
basin responses through physically distinct amplitude-modulation resonances determined by each
basin's geographic geometry and bathymetry.

---

## Cross-Validation Design: The 1940–1970 Test Window

### Why 1940–1970?

The interval 1940–1970 is scientifically important for several reasons:

1. **Data density**: The PSMSL (Permanent Service for Mean Sea Level) records are essentially
   complete for European stations from the 1930s onward.  Most long Baltic and North Sea
   records extend back to the 1860s–1890s, giving the model at least 50–80 years of training
   data before the holdout begins.

2. **Climatic interest**: The 30-year window spans the latter part of a documented global
   mean sea level acceleration that began ~1920, the wartime anomalies of 1940–1945, and
   the beginning of the modern satellite-era precursor trends.  It also straddles the 1940s
   warm AMO phase and the subsequent mid-century cooling.

3. **Independence guarantee**: Models trained before 1940 have no access to post-WWII
   instrumental improvements, data recovery campaigns, or station relocations—the holdout
   data are genuinely out-of-sample.

### Reading the Results Files

Each subdirectory contains:
- **`lte_results.csv`** — columns: decimal year, model output (mm), observed anomaly (mm),
  manifold amplitude, [other diagnostics].  The **Pearson correlation between columns 2 and 3
  restricted to 1940.0–1970.0** is the primary validation metric used throughout this
  narrative.
- **`<site>site1940-1970.png`** — time-series overlay of the model (blue/orange) against the
  observed record (black/grey) for the 1940–1970 window.

---

## Results Overview

The table below summarises all 79 experiments ordered by 1940–1970 validation correlation.
The "Full-record r" column gives the correlation over the entire training+test span for
context.

| Site | Location | r (1940–70) | r (full) | Record span |
|------|----------|-------------|----------|-------------|
| 11   | Warnemünde 2, Germany | **0.721** | 0.823 | 1855–2024 |
| 88   | Ratan, Sweden | **0.707** | 0.727 | 1892–2025 |
| 68   | Landsort, Sweden | **0.701** | 0.745 | 1887–2007 |
| 78   | Stockholm, Sweden | **0.698** | 0.751 | 1889–2025 |
| nino4 | NINO4 (W. Pacific) | **0.639** | 0.803 | 1880–2024 |
| 72   | Ystad, Sweden | **0.638** | 0.726 | 1887–1982 |
| 203  | Furuögrund, Sweden | **0.630** | 0.736 | 1916–2025 |
| 73   | Varberg, Sweden | **0.627** | 0.734 | 1887–1982 |
| 179  | Smögen, Sweden | **0.626** | 0.742 | 1911–2025 |
| 7    | Cuxhaven 2, Germany | 0.619 | 0.685 | 1843–2024 |
| 82   | Copenhagen, Denmark | 0.618 | 0.739 | 1889–2018 |
| ic3tsfc | ic3tsfc (recon. Pacific SST) | 0.616 | 0.780 | 1880–2016 |
| 14   | Helsinki, Finland | 0.610 | 0.708 | 1879–2025 |
| 70   | Kungsholmsfort, Sweden | 0.608 | 0.735 | 1887–2025 |
| 330  | Klagshamn, Sweden | 0.602 | 0.762 | 1930–2025 |
| 81   | Fredericia, Denmark | 0.601 | 0.740 | 1911–2025 |
| 172  | Mäntyluoto, Finland | 0.598 | 0.707 | 1911–2025 |
| 71   | Hanko/Hangö, Finland | 0.574 | 0.683 | 1888–2025 |
| 69   | Öland N. Cape, Sweden | 0.568 | 0.735 | 1887–2025 |
| 285  | Kaskinen/Kasko, Finland | 0.556 | 0.705 | 1927–2025 |
| 229  | Kemi, Finland | 0.548 | 0.695 | 1920–2024 |
| 57   | Vaasa/Vasa, Finland | 0.547 | 0.692 | 1884–2023 |
| amo  | AMO (N. Atlantic) | 0.543 | 0.791 | 1880–2023 |
| 95   | North Shields, UK | 0.536 | 0.672 | 1895–2025 |
| 98   | Slipshavn, Denmark | 0.531 | 0.733 | 1896–2018 |
| 8    | Wismar 2, Germany | 0.530 | 0.748 | 1849–2024 |
| 240  | Raahe/Brahestad, Finland | 0.528 | 0.667 | 1923–2025 |
| 79   | Oulu/Uleåborg, Finland | 0.526 | 0.706 | 1889–2025 |
| 24   | Delfzijl, Netherlands | 0.518 | 0.679 | 1865–2024 |
| 113  | Korsør, Denmark | 0.509 | 0.697 | 1897–2018 |
| nao  | NAO | 0.501 | 0.775 | 1880–2023 |
| 80   | Esbjerg, Denmark | 0.500 | 0.694 | 1889–2018 |
| 256  | La Jolla, CA, USA | 0.465 | 0.680 | 1925–2025 |
| nino34 | NINO3.4 | 0.449 | 0.698 | 1880–2024 |
| 2    | Świnoujście, Poland | 0.446 | 0.641 | 1811–2000 |
| 89   | Hirtshals, Denmark | 0.431 | 0.679 | 1892–2018 |
| 91   | Frederikshavn, Denmark | 0.430 | 0.684 | 1894–2018 |
| 62   | Oslo, Norway | 0.429 | 0.642 | 1886–2025 |
| 42   | Sevastopol, Ukraine | 0.413 | 0.698 | 1910–1995 |
| 25   | Harlingen, Netherlands | 0.413 | 0.681 | 1865–2024 |
| denison | Ft Denison, Australia | 0.411 | 0.540 | 1888–2024 |
| 225  | Ketchikan, AK, USA | 0.409 | 0.623 | 1919–2025 |
| emi  | EMI (ENSO Modoki) | 0.402 | 0.665 | 1880–2019 |
| 236  | West-Terschelling, Neth. | 0.401 | 0.641 | 1921–2024 |
| iodw | IOD West pole | 0.398 | 0.882 | 1880–2022 |
| 23   | Den Helder, Netherlands | 0.398 | 0.680 | 1865–2024 |
| 47   | Stavanger, Norway | 0.387 | 0.654 | 1919–2025 |
| nino3 | NINO3 | 0.382 | 0.701 | 1880–2024 |
| 1    | Brest, France | 0.377 | 0.679 | 1807–2025 |
| 58   | Bergen, Norway | 0.373 | 0.686 | 1915–2025 |
| 155  | Honolulu, HI, USA | 0.372 | 0.678 | 1905–2024 |
| 20   | Vlissingen, Netherlands | 0.371 | 0.632 | 1862–2024 |
| 111  | Fremantle, Australia | 0.367 | 0.688 | 1897–2023 |
| iode | IOD East pole | 0.364 | 0.782 | 1880–2022 |
| 161  | Galveston II, TX, USA | 0.363 | 0.667 | 1908–2025 |
| pdo  | PDO | 0.362 | 0.655 | 1880–2022 |
| 76   | Aarhus, Denmark | 0.348 | 0.654 | 1889–2018 |
| m4   | M4 (NPGO) | 0.341 | 0.715 | 1910–2024 |
| 154  | Trieste, Italy | 0.331 | 0.517 | 1877–2021 |
| 32   | IJmuiden, Netherlands | 0.309 | 0.617 | 1871–2024 |
| darwin1880 | Darwin SOI | 0.308 | 0.648 | 1880–2024 |
| 41   | Poti, Georgia | 0.306 | 0.651 | 1874–2020 |
| 239  | Turku/Åbo, Finland | 0.306 | 0.617 | 1922–2024 |
| 5    | Holyhead, UK | 0.303 | 0.297 | 1938–2025 |
| 246  | Pensacola, FL, USA | 0.268 | 0.663 | 1923–2025 |
| 127  | Seattle, WA, USA | 0.233 | 0.570 | 1899–2025 |
| 22   | Hoek van Holland, Neth. | 0.222 | 0.610 | 1864–2024 |
| brestexcl | Brest from 1880 | 0.220 | 0.664 | 1880–2020 |
| nino12 | NINO1+2 (E. Pacific) | 0.218 | 0.609 | 1880–2024 |
| 245  | Los Angeles, CA, USA | 0.214 | 0.697 | 1924–2025 |
| 183  | Portland, ME, USA | 0.212 | 0.661 | 1912–2025 |
| m6   | M6 (Atlantic Niño) | 0.212 | 0.582 | 1910–2024 |
| 302  | Tregde, Norway | 0.208 | 0.640 | 1928–2025 |
| 249  | Föglö/Degerby, Åland Is. | 0.201 | 0.577 | 1924–2025 |
| 202  | Newlyn, UK | 0.175 | 0.466 | 1915–2025 |
| tsa  | TSA (Trop. S. Atlantic) | 0.166 | 0.772 | 1880–2025 |
| iod  | IOD (composite) | 0.156 | 0.669 | 1880–2022 |
| tna  | TNA (Trop. N. Atlantic) | 0.128 | 0.761 | 1880–2025 |
| 234  | Charleston I, SC, USA | 0.119 | 0.559 | 1922–2024 |

---

## Baltic Sea and Gulf of Bothnia: The Highest-Performing Cluster

The single most striking pattern in the entire validation study is the **dominance of Baltic
Sea and Gulf of Bothnia tide gauges** at the top of the correlation table.  Sites 11, 88, 68,
78, 72, 203, 73, 179, 82, 14, 70, 330, 81, 172, 71, 69, 285, 229, 57, 240, 79, 8, 113, and
98—all located around the nearly enclosed Baltic basin—achieve 1940–1970 validation
correlations between **0.50 and 0.72**, uniformly higher than any site outside the region
except two Pacific climate indices (NINO4 and ic3tsfc).

### Warnemünde 2 (Site 11): The Reference Station

Warnemünde achieves the best validation score of any tide gauge in the study (r = 0.721).
With a record stretching from 1855 to 2024, it provides 85 years of training data before the
holdout begins—enough for the LTE solver to lock onto the dominant basin wavenumbers with
high confidence.  The site's LTE modulation settles on ltep₁ = 23.24 (the dominant
low-wavenumber Baltic basin resonance) and ltep₂ = 9.47 (the secondary basin mode), with
`harm 8` and `harm 24` selecting specific LTE octave harmonics.  The full-record correlation
(r = 0.823) is the highest achieved by any tide gauge in the study, confirming that
Warnemünde's sea-level signal is unusually coherent with the manifold across all amplitude
levels.

### Stockholm (Site 78) and Landsort (Site 68): The Oldest Continuous Baltic Records

Stockholm (1889–2025, r_val = 0.698) and Landsort (1887–2007, r_val = 0.701) are paired
siblings that together span more than 135 years.  Landsort hosts Sweden's oldest operating
lighthouse (1669) and its tide gauge is one of the PSMSL's reference quality stations.  Both
converge on very similar ltep parameters (ltep₁ = 23.0–23.3, ltep₂ = 9.47), strongly
suggesting that the Baltic tidal manifold amplitude response is a property of the enclosed
sea as a whole rather than of individual harbours.  The slightly higher validation score at
Landsort compared with Stockholm may reflect that Landsort is farther from the urban
heat/runoff effects of the capital.

### Ratan (Site 88) and the Gulf of Bothnia Gradient

Ratan (63.99°N, r_val = 0.707) is the northernmost site achieving super-0.70 validation.  It
sits at a critical geographic transition zone: to its north, postglacial land uplift (GIA)
accelerates and begins to dominate the sea-level signal; to its south, the Baltic proper's
tidal-meteorological forcing prevails.  Its remarkably high validation score—equal to that of
the Stockholm cluster—implies that the GIA-corrected residual signal in the northern Bothnian
Gulf still preserves strong lunisolar amplitude imprinting.  Moving northward to Furuögrund
(203, 64.92°N, r_val = 0.630), Raahe (240, 64.67°N, r_val = 0.528), Oulu (79, 65.04°N,
r_val = 0.526), Kemi (229, 65.67°N, r_val = 0.548), and ultimately Kaskinen (285, 62.34°N,
r_val = 0.556), the validation correlations form a physically coherent gradient that reflects
the increasing importance of GIA and freshwater fluxes at higher latitudes.

### Finnish Gulf of Finland Stations

Helsinki (14), Hanko/Hangö (71), and the Åland cluster (239 Turku, 249 Föglö) straddle the
entrance to the Gulf of Finland and the Archipelago Sea.  Helsinki (r_val = 0.610) performs
comparably to the Swedish Baltic stations, sharing almost the same ltep₁ ≈ −23.3 (the
negative sign indicates that the LTE modulation waveform is in anti-phase relative to the
Swedish convention—physically equivalent but with opposite amplitude polarity).  Hanko
(r_val = 0.574) is slightly lower, perhaps reflecting the bay's shallower geometry.  Turku
(r_val = 0.306) and Föglö (r_val = 0.201) are conspicuously weaker; both records start in
the early 1920s, leaving only about 18–20 years of training data before the holdout.  This
is likely a **training-length effect**: the solver does not have enough pre-1940 data to
stably constrain both wavenumber parameters simultaneously.

### Connecting the Baltic Cluster

The near-universal convergence of Baltic stations on ltep₁ ≈ 23.2 and ltep₂ ≈ 9.47 is
the most important structural finding of this study.  These two amplitude-modulation
wavenumber parameters are interpreted as the two lowest resonant modes of the Baltic–North
Sea enclosed basin: the ltep ≈ 23 mode as the fundamental long-scale standing-wave resonance
of the basin geometry, and the ltep ≈ 9.5 mode as its first spatial harmonic.  Together they
impose a characteristic multi-decadal beating pattern on the manifold amplitude that
reproduces the observed multi-decadal oscillations of Baltic sea level documented in the
tide-gauge literature.  The fact that the optimiser independently selects nearly identical
ltep values at sites separated by more than 1,000 km—from Klagshamn (55.5°N) to Kemi
(65.7°N)—is strong evidence that the recovered modulation is a property of the enclosed
Baltic basin as a whole rather than an overfitting artefact at individual sites.

---

## North Sea, Skagerrak, and Kattegat

### German Bight: Cuxhaven and Wismar

Cuxhaven (7, r_val = 0.619) and Wismar (8, r_val = 0.530) are North Sea / western Baltic
border sites.  Cuxhaven, situated at the mouth of the Elbe (53.87°N, 8.72°E), records one of
the longest uninterrupted European tide gauge series (from 1843).  Its validation performance
matches the lower end of the Baltic cluster even though it is exposed to North Sea tidal
dynamics—evidence that the Baltic basin amplitude resonance (ltep₁ ≈ 23) extends through the
German Bight.  Wismar (53.9°N, 11.5°E), on the southern Baltic coast just inside the Fehmarn
Belt, is geographically near Warnemünde but achieves a somewhat lower validation score,
consistent with its more complex local topography at the Baltic–North Sea transition.

### Danish Straits: Copenhagen, Fredericia, Korsør, Slipshavn

The Danish inner waters—spanning the Little Belt, Great Belt, and Øresund—form a geographic
transition zone between North Sea and Baltic proper.  Copenhagen (82, r_val = 0.618),
Fredericia (81, r_val = 0.601), Korsør (113, r_val = 0.509), and Slipshavn (98, r_val = 0.531)
all perform well above the North Sea average.  This consistency across four topographically
distinct straits implies that the shared ltep ≈ 23 Baltic basin amplitude mode is already
established at the point where Atlantic water enters the Baltic through the Danish transition
zone.

### Kattegat and Skagerrak Margins

Moving from the Kattegat into the Skagerrak, validation performance drops noticeably.
Aarhus (76, r_val = 0.348) and Esbjerg (80, r_val = 0.500) on Jutland's Kattegat/North Sea
coasts bracket this transition.  Esbjerg performs better than Aarhus despite facing directly
into the North Sea; its harbour geometry may filter Atlantic surge variability.  Hirtshals (89,
r_val = 0.431) and Frederikshavn (91, r_val = 0.430) at the Skagerrak–Kattegat junction show
a moderate but reproducible step down in validation skill relative to the central Danish Belt
sites, consistent with growing exposure to open-ocean meteorological forcing not captured by
the tidal manifold alone.

---

## Norwegian Atlantic Coast

### Bergen and Stavanger

Bergen (58, r_val = 0.373) and Stavanger (47, r_val = 0.387) are the primary Norwegian
Atlantic coast stations with records beginning in 1915 and 1919 respectively—both short enough
to limit pre-holdout training data to roughly 25 years.  Despite their short training period
both achieve validation correlations above 0.37; had longer records been available, performance
likely would have been higher.  Bergen's ltep₁ = 23.12 is indistinguishable from the Baltic
template, suggesting that the North Sea shelf acts as a resonant geographic waveguide
channelling the same basin amplitude mode up the Norwegian coast.

### Oslo (Site 62)

Oslo (r_val = 0.429) sits at the head of the Oslofjord (59.9°N, 10.7°E), a long narrow fjord
that acts as a selective geographic filter for sea-level variability.  The Oslofjord
integrates North Sea signals over a several-hundred-km fetch before they register at the
gauge, which may explain why Oslo's performance falls slightly below the Skagerrak sites
despite its record extending back to 1885.

### Smögen and North Shields: Cross-Skagerrak Coherence

Smögen (179, r_val = 0.626) on Sweden's Skagerrak coast is a standout: it achieves
essentially the same validation correlation as the top Baltic cluster, with ltep₁ = 22.88 and
ltep₂ = 9.45 consistent with the Baltic template.  North Shields (95, r_val = 0.536) on
England's North Sea coast (55.0°N, 1.4°W) similarly carries ltep₁ = 23.04 and ltep₂ = 9.41.
The strong performance of these two open North Sea sites—separated by the full width of the
Skagerrak—demonstrates that the ltep ≈ 23 basin amplitude mode is not a Baltic-specific
phenomenon but is detectable across the wider European shelf, consistent with the geographic
extent of the Baltic–North Sea enclosed basin system.

### Tregde (Site 302)

Tregde (r_val = 0.208), a sheltered harbour on Norway's Skagerrak coast (58.0°N, 7.6°E), is
one of the weaker Norwegian performers.  The record begins only in 1927.8, leaving fewer than
13 years of pre-holdout training data.  Its full-record correlation (0.640) is healthy,
indicating that the manifold fit is physically sensible over the full period; the weak
validation score reflects the underconstrained pre-1940 calibration more than any fundamental
mismatch.

---

## UK and Irish Sea

### Newlyn (Site 202): Britain's Sea-Level Datum

Newlyn (r_val = 0.175) is the reference station for Ordnance Datum Newlyn, the benchmark for
all elevations in Great Britain, yet it achieves the *lowest* validation correlation of any
European site in the study.  This is a striking and potentially important result.  The Newlyn
record is of the highest quality (PSMSL quality-A, very few gaps), so the low score is not an
artefact.  Its ltep₁ = −10.22 and ltep₂ = 23.59 are atypical: the primary resonant mode is
a negative ~10 wavenumber parameter rather than the ltep ≈ 23 template shared by all
neighbouring sites.  The Newlyn tide gauge (50.1°N, 5.5°W) is exposed to the full fetch of
the north Atlantic via the English Channel—a geography that allows storm surges, long-period
seiches, and Atlantic meridional overturning variability to contribute amplitude variance in
ways that the tidal manifold cannot replicate.

### Holyhead (Site 5): Irish Sea

Holyhead (r_val = 0.303) achieves a full-record correlation of only 0.297—the only site in
the study where the validation correlation *exceeds* the full-record correlation.  Its record
starts in 1938, meaning that only two years of data precede the holdout; the solver is
effectively operating in pure prediction mode.  The positive (if modest) validation
correlation under these extreme conditions is a compelling demonstration that the manifold
carries genuine predictive information even without a long pre-holdout calibration window.

---

## North Atlantic: French, German, and Dutch Open Coast

### Brest (Site 1 and brestexcl)

Brest (1, r_val = 0.377) holds the world's longest continuous tide gauge record, stretching
from 1807 to the present.  With 133 years of pre-holdout training data, the model is as well
calibrated as physically possible, yet the validation correlation is only moderate.  The
Brest (from-1880) version (`brestexcl`, r_val = 0.220) performs substantially worse with its
shorter pre-holdout window, confirming the training-length dependence already noted for the
Baltic.

The fact that Brest's enormous training dataset yields only r_val = 0.377 has a clear
geographic interpretation: Brest (48.4°N, 4.5°W) is exposed to the full North Atlantic,
where Atlantic Multidecadal Oscillation (AMO) variability, Gulf Stream meanders, and
storm-track changes contribute amplitude variance on timescales that are orthogonal to the
lunisolar forcing.  These climate-driven contributions are geographically distinct from the
enclosed-basin resonance mechanism and no amount of training data will improve their fit.

### Dutch Wadden Sea Cluster

The five principal Dutch gauges—Vlissingen (20, r_val = 0.371), Hoek van Holland (22,
r_val = 0.222), Den Helder (23, r_val = 0.398), Delfzijl (24, r_val = 0.518), West-
Terschelling (236, r_val = 0.401), and Harlingen (25, r_val = 0.413)—span a wide performance
range.  Delfzijl (53.3°N, 7.0°E), at the mouth of the Ems estuary adjacent to the Groningen
gas fields, is the best performer; its ltep parameters match the Baltic template, suggesting
that at the eastern edge of the Wadden Sea, the Baltic-type basin amplitude mode (routed
around Jutland via the Skagerrak/Kattegat) already dominates over the open North Sea signal.
Hoek van Holland (22, r_val = 0.222) is unexpectedly weak; with its ltep₁ = −22.37 (a
sign-flipped mode), it may be capturing a North Sea channel resonance that is in anti-phase
with the Baltic convention—a consequence of its geographic position at the Rhine outflow.

IJmuiden (32, r_val = 0.309) anchors the North Sea Canal and has the second-lowest Pearson r
of the Dutch cluster.  Its ltep₁ = 23.84 is the *longest* primary wavenumber of any Dutch
station, shifted by nearly 0.7 units from Delfzijl's 23.12—a small but apparently
consequential geographic offset in amplitude response for the 1940–1970 window.

---

## Black Sea and Adriatic

### Sevastopol (Site 42) and Poti (Site 41)

Sevastopol (r_val = 0.413) and Poti (r_val = 0.306) are the only Black Sea representatives
in the study.  Both achieve moderate validation despite the Black Sea being effectively an
enclosed sea with minimal direct tidal forcing—exactly the kind of geography where the common
tidal manifold might be expected to fail.  Their success suggests that even in the absence of
direct tidal forcing, the lunisolar-driven atmospheric circulation modulated at the
Baltic-equivalent amplitude scale (ltep₁ = 23.23 at Sevastopol) still impresses a coherent
signal on Black Sea sea level.  Sevastopol's ltep₁ = 23.23 is indistinguishable from the
Baltic template in numeric value, suggesting that the Black Sea's enclosed geographic geometry
supports a similar fundamental resonant mode.  Poti's record begins in 1874 at the eastern
corner of the Black Sea (42.1°N, 41.7°E) where the Rioni River outflow introduces substantial
freshwater variability, which likely explains its lower validation performance compared with
Sevastopol.

### Trieste (Site 154): The Adriatic

Trieste (r_val = 0.331) at the head of the Adriatic Sea (45.6°N, 13.8°E) is a moderately
performing site.  Its ltep₁ = −20.73 is notably different from the ~23 template: the Adriatic
is a resonant elongated channel with its own geographic natural frequencies, and the LTE
amplitude modulation appears to have correctly identified this regional departure from the
pan-European Baltic pattern.  The negative sign and the distinct wavenumber value of 20.7
reflect the Adriatic's different enclosed-basin geometry and depth structure compared with
the Baltic–North Sea system.

---

## Atlantic Coast of North America

### New England: Portland, Maine (Site 183)

Portland, ME (r_val = 0.212) spans 1912–2025 but is effectively an open-ocean site (43.7°N,
70.2°W) exposed to Gulf Stream dynamics and the North Atlantic Oscillation.  Its full-record
correlation (0.661) confirms that the model captures the long-term amplitude envelope well;
it is specifically the 1940–1970 period that is difficult.  This is the geographic window
where Gulf Stream variability associated with the mid-century AMO transition—a period of
enhanced decadal climate noise—dominates over the tidal forcing.

### Mid-Atlantic and Gulf of Mexico: Charleston (Site 234), Galveston (Site 161),
### Pensacola (Site 246)

Charleston (r_val = 0.119) is the weakest-performing tide gauge in the study.  The
1940–1970 period at Charleston (32.8°N, 79.9°W) is complicated by the acceleration of
subsidence from groundwater extraction, which was ongoing through mid-century and adds a
non-stationary secular component that the linear background term (`bg`) cannot fully absorb.
Galveston (161, r_val = 0.363) and Pensacola (246, r_val = 0.268) similarly suffer from
subsidence, though Pensacola's longer pre-holdout window (1923 start) provides somewhat
better calibration.  The full-record correlations for both Gulf sites (0.663–0.667) are
substantially higher than the holdout, pointing again to a geographic-regime shift during
1940–1970 rather than a fundamental model failure.

---

## Pacific Rim

### US West Coast: La Jolla, Los Angeles, Seattle

La Jolla (256, r_val = 0.465), Los Angeles (245, r_val = 0.214), and Seattle (127,
r_val = 0.233) show a wide spread despite occupying similar North Pacific latitudes.  La Jolla
(32.9°N, 117.3°W) achieves the highest validation correlation of the three with ltep₁ = 7.08—
a sub-decadal wavenumber absent from Atlantic sites, likely reflecting the ENSO-related
amplitude modulation of the eastern Pacific at this basin-shelf geometry.  Los Angeles
(ltep₁ = −4.10) carries a smaller wavenumber that overlaps with the ENSO biennial amplitude
scale, while Seattle (ltep₁ = 24.11) converges on a value close to the Baltic template but
still achieves only moderate validation.  Seattle's poor 1940–1970 performance is consistent
with documented PDO-forced anomalies on the Washington/British Columbia coast during this
geographic period.

### Ketchikan, Alaska (Site 225)

Ketchikan (r_val = 0.409, 55.3°N, 131.6°W) performs better than any other Pacific US station
in the 1940–1970 validation.  Its record begins in 1919, providing 21 years of pre-holdout
training.  Ketchikan's ltep₁ = 23.40 aligns with the pan-European Baltic template, suggesting
that the same fundamental geographic basin resonance mode active in the Baltic–North Sea
system is also present in the North Pacific continental shelf geometry at this latitude.

### Honolulu, Hawaii (Site 155)

Honolulu (r_val = 0.372, 21.3°N, 157.9°W) is a mid-Pacific open-ocean site whose sea level
is closely tied to ENSO and the PDO.  Its ltep₁ = 22.32 is close to but slightly shorter
than the European ltep ≈ 23 template, and its secondary mode ltep₂ = 9.46 matches the Baltic
basin second harmonic.  The moderate validation performance indicates that the tidal manifold
captures a meaningful fraction of the Pacific interannual-to-decadal amplitude variance even
at this geographically isolated open-ocean site.

**Caveat and area for further study**: The ltep₁ ≈ 22.3, ltep₂ ≈ 9.46 parameter combination
recovered for Honolulu is essentially the Baltic template, which is unexpected for a tropical
mid-Pacific island site that is geographically remote from the Baltic–North Sea enclosed basin
system.  Based on Hawaii's known oceanographic setting—sea level at Honolulu is strongly
modulated by ENSO-related warm-pool variability advected from the western equatorial Pacific—
one would *a priori* anticipate a **NINO4-like wavenumber structure** (ltep₁ ≈ 12, the
resonant mode of the western Pacific warm-pool basin) rather than a Baltic-like mode near 23.
The fact that the optimiser instead selects ltep ≈ 23 suggests either that (a) the short
training window (1905–1939, only 35 years) is insufficient to distinguish the NINO4 and Baltic
modes whose multi-decadal beating patterns can appear similar on short records, or (b) the
secular background component shared with the Baltic template is absorbing a disproportionate
share of Hawaii's sea-level variance, masking the true ENSO-amplitude wavenumber structure.
This site is therefore flagged as a **priority for further investigation**: re-running the
optimisation with the ltep parameters explicitly constrained to the NINO4 basin range
(ltep ≈ 10–14) and with a detrended training record would clarify whether Hawaii genuinely
responds at the Baltic basin wavenumber or whether the current fit is a consequence of secular
confounding.

### Fort Denison, Sydney (denison)

Fort Denison (r_val = 0.411) in Sydney Harbour (33.9°S, 151.2°E) is the only Southern
Hemisphere MSL station in the study.  Its full-record correlation (0.540) is lower than most
European sites, consistent with the Southern Ocean's more chaotic mid-latitude circulation.
The `ltep` structure for Fort Denison is unusual: ltep₁ = 0.01 ≈ 0 (effectively vanishing
first mode) and ltep₂ = 23.71—a configuration that shifts most of the LTE amplitude
modulation into the second wavenumber parameter.  The positive validation r (0.411) suggests
that the ltep ≈ 23 basin mode is still identifiable in the Southern Hemisphere record even
at this geographic distance from the Baltic system.

**Fremantle-like NINO4 character**: Despite its moderate cross-validation score, Fort Denison
actually exhibits a wavenumber structure that can be interpreted as possessing a **Fremantle /
NINO4 geographic character**.  Fremantle (ltep₁ = 12.02) is located on Australia's Indian
Ocean coast and is directly connected to the central Pacific ENSO signal via the Leeuwin
Current.  Fort Denison on the Pacific coast of Australia (eastern seaboard, 33.9°S) receives
a related but Pacific-routed ENSO signal through the East Australian Current.  The vanishing
ltep₁ at Fort Denison, combined with the secondary ltep₂ = 23.71 that is closest to the
Baltic template, may reflect the optimiser's difficulty in cleanly separating the NINO4-type
equatorial Pacific amplitude mode (ltep ≈ 12) from the background secular drift when only
~52 years of pre-holdout data are available (1888–1939).  The moderate r_val = 0.411 is
therefore more encouraging than it first appears: it indicates that the manifold is genuinely
capturing ENSO-linked amplitude variability at Fort Denison consistent with an eastern Pacific
boundary pathway, analogous to Fremantle on the western side.  **This site is a candidate for
further investigation** with the ltep search constrained to the NINO4 basin range (ltep ≈ 12)
to determine whether an explicit ENSO-amplitude fit would improve or maintain the
cross-validation performance observed here.

### Fremantle, Australia (Site 111)

Fremantle (r_val = 0.367, 32.1°S, 115.7°E) on the Indian Ocean coast of Western Australia
is strongly influenced by the Leeuwin Current, an anomalously warm southward-flowing
boundary current that carries ENSO-related signals from the tropical Pacific into the
Australian Bight.  Its ltep₁ = 12.02 matches the NINO4 and NINO3.4 ltep₁ values (12.1 and
12.9 respectively), providing a direct **geographic connection** between Fremantle's
sea-level amplitude variability and the central Pacific ENSO signal: the Indo-Pacific basin
geometry at these scales supports a resonant amplitude modulation near ltep ≈ 12 that is
identifiable both in the Leeuwin Current boundary signal at Fremantle and in the Pacific
ENSO indices.

---

## Climate Indices

### NINO4 and ic3tsfc: The Top-Performing Non-Tide-Gauge Experiments

NINO4 (r_val = 0.639) is the best-performing climate index in the study, surpassed only by
the top four Baltic tide gauges.  The NINO4 region (equatorial central Pacific, 5°N–5°S,
160°E–150°W) represents the western warm pool of the Pacific El Niño system.  Its
ltep₁ = 12.13 is the dominant resonant wavenumber of the western equatorial Pacific basin
geometry, reflecting the amplitude-level response of the warm-pool basin to the lunisolar
manifold at this geographic scale.  The remarkably high full-record correlation (r = 0.803)
drops by only 0.16 in the holdout (r_val = 0.639), indicating excellent model stability.

The ic3tsfc reconstruction (r_val = 0.616) is an instrumental-era SST reconstruction for the
tropical Pacific with ltep₁ = −12.22—nearly identical in amplitude-scale to NINO4's 12.13
but with flipped sign.  The sign flip corresponds to an anti-phase amplitude modulation
waveform, essentially selecting the inverse ENSO mode that peaks in central Pacific cooling
(La Niña–type) when NINO4 peaks in warming.  The close agreement between NINO4 and ic3tsfc
in both validation correlation and LTE wavenumber structure provides internal consistency
evidence that the model is resolving real geographic amplitude variability.

### AMO (amo): Multidecadal Atlantic Memory

The Atlantic Multidecadal Oscillation (r_val = 0.543) achieves the third-highest validation
correlation among climate indices.  Its unusual ltep structure—ltep₁ = −7.17 and
ltep₂ = −7.22, both negative and nearly equal—produces a strongly symmetric amplitude
modulation waveform, consistent with the AMO's known quasi-decadal pulsing within its
multidecadal envelope.  The full-record AMO correlation (0.791) drops by about 0.25 in the
holdout, a gap that likely reflects the AMO's known regime shifts at ~1965 and ~1995 that
are imperfectly captured by a stationary LTE amplitude modulation.

A key geographic connection: the AMO index encapsulates much of the mechanism that limits
Atlantic tide gauge validation.  Sites like Brest (1, r_val = 0.377), Newlyn (202, r_val = 0.175),
and Charleston (234, r_val = 0.119) all occupy basins where AMO amplitude variance dominates
on the same scales at which the tidal manifold also operates.  The AMO itself achieves
moderate validation (0.543) but does so through a *different* LTE wavenumber parameter than
the tide gauges—ltep ≈ 7 rather than ltep ≈ 23—suggesting that the AMO's North Atlantic
basin geometry and the European tide gauges' Baltic basin geometry are related but distinct
geographic systems.

### NAO (nao): The Atmospheric Bridge to Baltic Sea Level

The North Atlantic Oscillation (r_val = 0.501) is the dominant atmospheric mode driving
European weather and sea level.  Its ltep₁ = 24.62 is the longest primary wavenumber of any
well-validated index in the study, slightly larger than the Baltic cluster ltep₁ ≈ 23.2.
The NAO's comparatively high validation score—despite being a pressure index rather than an
ocean variable—is consistent with its close **geographic relationship** to the Baltic Sea
region: Baltic MSL variability is strongly wind-driven, and the NAO centres of action
(Iceland Low / Azores High) directly modulate westerly wind stress over the North Sea and
the Baltic inflow through the Danish Straits.  The near-coincident ltep values of the NAO
(≈ 24.6) and the Baltic tide gauges (≈ 23.2) reflect that both are responding to the same
amplitude modulation of the lunisolar manifold as projected through the geographic basin
geometry of the Baltic–North Sea system.  The NAO is therefore geographically most naturally
connected to Baltic MSL sites rather than to any open-ocean Atlantic tide gauge.

### ENSO Structure: The Pacific Geographic Gradient

The four NINO indices reveal a systematic east-to-west geographic gradient in validation
performance:

| Index | Region | ltep₁ | r_val | r_full |
|-------|--------|--------|-------|--------|
| NINO4 | W. Pacific warm pool | 12.13 | 0.639 | 0.803 |
| NINO3.4 | Central-E. Pacific | 12.87 | 0.449 | 0.698 |
| NINO3 | E. Pacific | 94.3 | 0.382 | 0.701 |
| NINO1+2 | Far-E. Pacific | 94.0 | 0.218 | 0.609 |

The ltep₁ transition from ≈ 12 in the west to ≈ 94 in the east is dramatic.  The
western warm pool (NINO4 geographic region) is dominated by the ltep ≈ 12 resonant
wavenumber of the equatorial Pacific warm-pool basin geometry; the eastern boundary (NINO3,
NINO1+2 geographic regions) is dominated by a very large-scale mode that is nearly secular
on the amplitude scale of the model.  This is consistent with oceanographic understanding:
the eastern Pacific SST is strongly influenced by coastal upwelling driven by trade wind
stress and the annual ITCZ cycle—a geographic setting where the basin amplitude response
integrates over much larger spatial scales—while the western warm pool integrates
multi-decadal heat content anomalies that respond at the basin-scale resonant wavenumbers
of the equatorial Pacific.

The ENSO Modoki Index (EMI, r_val = 0.402) represents a distinct flavour of El Niño
characterised by central Pacific heating rather than eastern-boundary heating.  Its
ltep₁ = 63.29 sits between the ≈ 12 western-pool mode and the ≈ 94 eastern-boundary
mode, consistent with its intermediate geographic character.

### Indian Ocean Dipole: A Cautionary Tale

The IOD composite index (r_val = 0.156) and its east (iode, r_val = 0.364) and west
(iodw, r_val = 0.398) poles illustrate how a strong full-record correlation can coexist with
weak holdout validation.  The IOD west pole achieves r_full = 0.882—the highest full-record
correlation in the entire study—yet drops to r_val = 0.398 in the holdout.  The IOD east pole
similarly falls from r_full = 0.782 to r_val = 0.364.

This collapse is physically meaningful.  The IOD indices have a known strong covariance with
ENSO on seasonal-to-interannual timescales, but this covariance is not stationary: the
IOD–ENSO geographic coupling strengthened after ~1970 and was different during the 1940–1970
period.  The LTE model, trained on the full record (which includes post-1970 data), learns an
LTE wavenumber parameter tuned to the strong post-1970 IOD–ENSO coupling, which then
*underperforms* in the 1940–1970 window when that coupling was weaker.  The IOD results
therefore serve as an important validation integrity check: the drop in holdout performance
for the IOD indices is not a model failure but a correct detection of non-stationarity in the
index itself.

**Secular-trend caveat for iode and iodw**: An additional interpretive complication arises
because both IOD poles (iode and iodw) carry a pronounced **secular upward trend** over their
1880–2022 record.  When the optimiser fits the manifold over the full record, the secular ramp
embedded in the manifold's background acceleration (`bg`) can absorb this trend very
efficiently, producing a high r_full largely driven by the shared long-term drift rather than
by the cyclic interannual-to-decadal variability that the LTE modulation is designed to
capture.  The holdout window 1940–1970 is a relatively stationary segment of the record where
the secular slope is compressed into a modest fraction of the total variance; the model's
over-reliance on matching the trend during training is then *exposed* by the cross-validation,
and the drop in r from r_full to r_val is partly a consequence of this secular-fit bias rather
than a pure non-stationarity effect.  In other words, **the "poorer" cross-validation score
for IOD indices partly reflects the training procedure rewarding the secular trend at the
expense of the cyclic component that the holdout actually tests**.

### Tropical Atlantic Indices: TNA, TSA, M6

The Tropical North Atlantic (TNA, r_val = 0.128), Tropical South Atlantic (TSA, r_val = 0.166),
and Atlantic Niño (M6, r_val = 0.212) all exhibit the same pattern as the IOD: high full-record
correlations (r_full = 0.761, 0.772, 0.582 respectively) collapsing in the 1940–1970 holdout.
These Atlantic tropical SST indices are known to be sensitive to mid-century volcanic aerosol
forcing (the 1963 Agung and 1968 Fernandina eruptions straddle the holdout) as well as to AMO
phase transitions, neither of which is represented in the tidal manifold.  The geographic
contrast between the tropical Atlantic indices (r_val ≈ 0.13–0.21) and the tropical Pacific
NINO indices (r_val ≈ 0.22–0.64) quantifies the degree to which the Pacific basin geometry is
more strongly coupled to lunisolar amplitude forcing than the Atlantic on these scales.

**Secular-trend caveat for tsa and tna**: Like the IOD poles, both the Tropical South Atlantic
(TSA) and Tropical North Atlantic (TNA) indices possess a **strong monotonic warming trend**
over their 1880–2025 records.  During optimisation on the full record, the manifold's secular
background ramp can fit this trend with high fidelity, yielding inflated r_full values that
reflect shared long-term drift rather than reproduced cycle-to-cycle variability.  The 1940–1970
holdout exposes this weakness: the cyclic component of TNA and TSA variability that actually
varies within that stationary window is small relative to the secular offset, and the model—
having devoted much of its fitting capacity to the multi-decade trend—lacks residual skill to
predict the interannual variations that dominate the holdout metric.  Readers should therefore
interpret the low r_val for tna (0.128) and tsa (0.166) as partly reflecting this
**secular-trend fitting bias** rather than purely indicating that lunisolar forcing has no
relationship to tropical Atlantic variability.  A fairer assessment would detrend both the
model and observations before computing the holdout correlation, or restrict the training to a
sub-period equally affected by the trend.

### PDO and NPGO

The Pacific Decadal Oscillation (PDO, r_val = 0.362) with ltep₁ = 4.72 and the North
Pacific Gyre Oscillation (M4/NPGO, r_val = 0.341) with ltep₁ = −7.08 both achieve
modest validation.  The PDO's ltep₁ ≈ 4.7 reflects the amplitude response of the North
Pacific gyre to the ENSO-related lunisolar forcing at sub-decadal basin scales.  The NPGO's
ltep₁ ≈ −7.08 matches the AMO's ltep₁ ≈ −7.2 in scale, suggesting a shared sub-basin
amplitude mode connecting North Pacific gyre circulation and North Atlantic at similar
geographic scales—an inter-ocean teleconnection expressed through their respective basin
geometries.

### QBO (qbo30): Stratospheric Oddity

The Quasi-Biennial Oscillation at 30 hPa stands apart structurally from all other experiments.
Its `bg` = 5.26 yr⁻¹ (compared with 1.4 × 10⁻⁴ yr⁻¹ for all other sites) and
`delA` = −0.149 (negative, compared with +4.5 for all others) signal that the LTE solver
has fit the QBO's fast oscillation with a completely different manifold geometry.  No
validation PNG is present for the QBO in the 1940-1970 window because the QBO instrumental
record starts only in 1953; the `lte_results.csv` covers a different window.

---

## Discussion: Cross-Site Coherence and the Common Manifold

### The ltep₁ ≈ 23 Baltic Basin Wavenumber as a Pan-European Signal

The single most important cross-site finding is the remarkable convergence of the primary LTE
wavenumber parameter onto **ltep₁ ≈ 23.2 ± 1** across all Baltic, North Sea, and
Skagerrak/Kattegat stations, as well as several Norwegian and UK sites.  Within the LTE
framework this is the **fundamental low-wavenumber resonance of the Baltic–North Sea enclosed
basin**—a spatial amplitude-modulation resonance rather than a temporal oscillation (see
*Mathematical Geoenergy*, Wiley/AGU, 2019; and the TeX source at https://github.com/pukpr/gem).

The physical basis is geographic and geometric: the bathymetry, coastline, and depth structure
of the Baltic–North Sea system define a preferred spatial scale for standing-wave amplitude
modulation of the tidal manifold.  This has nothing to do with a temporal cycle of any
particular duration; it is about how the enclosed geometry of the basin projects the
manifold's amplitude onto observable sea-level variance.  Different geographic basins produce
different resonant ltep values:

- The **Baltic–North Sea** geographic basin settles on ltep₁ ≈ 23.2 (fundamental resonance)
  and ltep₂ ≈ 9.5 (spatial harmonic), consistently across 25+ sites from Klagshamn in
  southern Sweden to Kemi on the northern Gulf of Bothnia.
- The **western equatorial Pacific warm pool** (NINO4 index; Fremantle via the Leeuwin
  Current) resolves near ltep₁ ≈ 12, reflecting the larger geographic basin dimensions
  of the Pacific warm pool.
- The **Adriatic Sea** finds ltep₁ ≈ −20.7, consistent with the Adriatic's own lower
  resonant wavenumber relative to the broader European shelf—a consequence of the Adriatic's
  distinct elongated geometry.
- The **eastern Pacific coastal boundary** (NINO3, NINO1+2) produces ltep₁ ≈ 94, indicating
  that upwelling-dominated coastal settings at the South American boundary respond at much
  larger spatial amplitude scales.

The LTE transformation selects the resonant wavenumber that best projects the manifold onto
the observed sea-level variance at each geographic site.  The optimiser independently
recovering ltep₁ ≈ 23 at 25+ geographically separated European sites—spanning more than
1,000 km from Klagshamn (55.5°N) to Kemi (65.7°N)—is compelling evidence that this is a
genuine resonance of the Baltic–North Sea basin geometry rather than an overfitting artefact.

### The ltep₂ ≈ 9.47 Second Spatial Harmonic

The second wavenumber parameter ltep₂ ≈ 9.47 appears as the dominant secondary mode for
virtually every Baltic and North Sea gauge.  Its ubiquity is interpreted as the first spatial
harmonic of the ltep₁ ≈ 23 fundamental: just as a resonant cavity supports overtone modes,
the Baltic–North Sea system supports a higher-wavenumber amplitude-modulation mode near
half the fundamental scale.  The slight deviation from an exact ratio (23/9.47 ≈ 2.43 rather
than 2.0) is consistent with the resonance being geometrically determined, leading to mode
shapes that are set by the basin's irregular coastline and varying bathymetry rather than by
simple harmonic arithmetic.

### Why Validation Degrades from Baltic to Open Ocean

The systematic gradient in validation performance—from r_val > 0.70 in the Baltic to
r_val < 0.20 for some Atlantic and tropical indices—is physically interpretable as the
**ratio of tidal-manifold-driven variance to total variance** in each time series:

1. **Baltic Sea (enclosed, geographically constrained)**: Regional sea level is dominated by
   wind stress, river runoff, and atmospheric pressure integrated over multi-year timescales.
   The enclosed basin geometry supports strong amplitude-modulation resonances at ltep ≈ 23
   and ltep ≈ 9.5, efficiently excited by the lunisolar manifold.  → High r_val.

2. **North Sea shelf (partially open, moderate tides)**: The same ltep ≈ 23 basin resonance
   is present but competes with NAO-driven storm surge variability and North Atlantic water
   mass exchange entering through the northern margin.  → Intermediate r_val.

3. **Open Atlantic coast** (Brest, Newlyn, Portland ME): Fully exposed to AMO, Gulf Stream,
   and storm track variability.  The signal-to-noise ratio for the basin-resonance amplitude
   contribution is low.  → Low r_val.

4. **Tropical Pacific indices** (NINO4, NINO3.4): The western Pacific warm pool geographic
   basin supports a resonant amplitude-modulation at ltep ≈ 12 that the manifold captures
   well.  → Moderate-to-good r_val.

5. **Tropical Atlantic and IOD indices**: These exhibit time-varying relationships with
   ENSO and with volcanic/aerosol forcing that breaks stationarity in the 1940–1970 window.
   → Low r_val, high r_full.

### Training Length vs. Physical Mechanism

A secondary factor modulating validation performance is training-record length.  Sites where
the record begins less than ~20 years before 1940—Holyhead (5, starts 1938), Föglö (249,
starts 1924), Turku (239, starts 1922), West-Terschelling (236, starts 1921), Tregde (302,
starts 1928)—all show markedly lower validation scores than their geographical neighbours
with longer records.  This indicates a practical threshold: the LTE solver needs roughly
**20–30 years of calibration data** to constrain the LTE wavenumber phase with sufficient
precision to extrapolate into the holdout window.  The Baltic stations that were instrumented
before 1880 (Warnemünde 1855, Stockholm 1889, Helsinki 1879, Landsort 1887, Cuxhaven 1843,
Wismar 1849, Brest 1807) all benefit from 50–100 years of pre-holdout calibration, which
translates into the consistently high validation correlations observed.

---

## Connections Between Sites: The Manifold as a Unifying Framework

One of the key promises of the common latent manifold approach is that it makes
*cross-site comparison* meaningful in a way that site-by-site empirical modelling does not.
Because every experiment in this directory is built from the same underlying lunisolar forcing
(the LPAP sinusoid table is essentially identical across sites), differences in LTE
modulation parameters directly reflect differences in the *geographic basin geometry* through
which each location responds to the forcing at the amplitude level.

The following geographic connections emerge from the parameter survey:

1. **Baltic ↔ North Shields ↔ Smögen**: All three converge on ltep₁ ≈ 23, ltep₂ ≈ 9.47,
   demonstrating that the European shelf's dominant basin resonance spans from the central
   North Sea (North Shields, 55.0°N, England; Smögen, 58.4°N, Sweden) all the way into the
   innermost Baltic.  This geographic continuum reflects a shared shelf-basin geometry rather
   than any coincidence of temporal cycles.

2. **Fremantle ↔ NINO4 ↔ ic3tsfc**: All three share ltep₁ ≈ 12, linking Indian Ocean–Pacific
   boundary current variability at Fremantle (32.1°S, Western Australia, where the Leeuwin
   Current carries tropical Pacific signals southward) to the western Pacific warm pool
   NINO4 region (5°N–5°S, 160°E–150°W).  The shared ltep value reflects a common geographic
   basin scale of the Indo-Pacific warm pool system at these latitudes.

3. **NAO ↔ Baltic MSL sites**: The NAO ltep₁ ≈ 24.6 is the closest of all climate-index
   wavenumber parameters to the Baltic cluster ltep₁ ≈ 23.2.  This geographic connection is
   direct: the NAO pressure centres of action (Iceland Low, Azores High) straddle the North
   Sea and control westerly wind stress over the Baltic inflow through the Danish Straits.
   The near-coincident ltep values reflect that both the atmospheric NAO and the Baltic tide
   gauges are driven by the same basin geometry—the Baltic–North Sea enclosed system—and its
   amplitude-level response to the lunisolar manifold.  The NAO is therefore geographically
   more naturally connected to Baltic MSL variability than to any open-ocean Atlantic site.

4. **AMO ↔ NPGO (M4)**: The AMO (ltep₁ = −7.2; North Atlantic) and the NPGO (ltep₁ = −7.08;
   North Pacific) share a similar wavenumber scale, suggesting a shared sub-basin resonance
   connecting the Atlantic multidecadal circulation to the North Pacific gyre.  Both are
   geographically removed from the Baltic cluster (ltep₁ ≈ 23), consistent with the Atlantic
   and Pacific basin geometries being distinct from the European shelf.

5. **IOD west ↔ IOD east ↔ EMI**: All share ltep₁ near 1.9 or 63—far from the European ~23
   cluster—reflecting the Indian Ocean's distinct geographic basin geometry and its strongly
   interannual (ENSO-modulated) amplitude structure.

6. **NINO3 ↔ NINO1+2**: Both converge on ltep₁ ≈ 94, a near-secular amplitude scale,
   quantifying how the eastern Pacific coastal boundary (from Baja California to northern
   Chile) integrates multi-decadal climate forcing over a very large geographic domain.

7. **Brest ↔ Brest-excl**: The ltep₁ values differ between the full Brest record
   (ltep₁ = +9.57) and the Brest-from-1880 version (ltep₁ = −9.88).  The sign flip and
   near-identical magnitude point to a phase ambiguity in the secondary basin mode when the
   18th-century data are excluded, consistent with the known data quality issues in the
   pre-1880 Brest record.

---

## Conclusion

The February 2026 cross-validation experiments demonstrate that the GEM-LTE common latent
manifold framework achieves robust out-of-sample predictive skill on the 1940–1970 holdout
interval for the majority of the 79 tested sites.  The most important conclusions are:

1. **The Baltic Sea–North Sea European shelf system is the strongest domain for manifold
   validation**, with 20+ sites achieving r_val > 0.50 and the top four (Warnemünde, Ratan,
   Landsort, Stockholm) exceeding r_val = 0.70.  The recovery of a consistent ltep₁ ≈ 23.2
   across all of these geographically dispersed stations is the core empirical result of the
   study, interpreted as the fundamental low-wavenumber resonant amplitude mode of the
   Baltic–North Sea enclosed basin geometry.

2. **The tropical central Pacific (NINO4, ic3tsfc) achieves validation skill comparable to
   the best tide gauges**, establishing a direct geographic link between the Pacific warm pool's
   amplitude variability and the ltep₁ ≈ 12 resonant mode of the equatorial Pacific basin.

3. **The NAO's wavenumber parameter (ltep₁ ≈ 24.6) is closest to the Baltic ltep₁ ≈ 23.2**,
   consistent with the NAO being geographically most directly coupled to Baltic MSL variability
   through its modulation of westerly wind stress over the North Sea and the Baltic inflow
   through the Danish Straits.  Of all climate indices tested, the NAO is the one most
   naturally linked to the Baltic tide gauge sites.

4. **Full-record versus holdout divergence is diagnostically informative**: indices where
   r_full >> r_val (IOD, tropical Atlantic) signal non-stationary coupling that invalidates
   stationarity assumptions, while sites where r_full ≈ r_val (most Baltic gauges) signal
   robust, physically stationary amplitude-modulation coupling between the basin geometry
   and the manifold.  A further mechanism producing large r_full/r_val gaps is **secular
   trend confounding**: indices with a strong monotonic trend (iode, iodw, tna, tsa) allow
   the manifold's background acceleration to absorb much of the full-record variance, inflating
   r_full at the expense of cyclic fit quality, which the detrended 1940–1970 holdout then
   exposes.

5. **Training-record length imposes a practical constraint of ~20–30 years** before the
   holdout for reliable LTE amplitude-phase calibration.

6. **The common manifold is not merely a fitting device** — its ability to reproduce the
   1940–1970 interval *without ever seeing that data* provides the strongest evidence yet that
   the ltep₁ ≈ 23 and ltep₂ ≈ 9.5 wavenumber parameters recovered by the model reflect
   genuine basin-geometry resonances rather than curve-fitting to the observed record.  The
   LTE modulation operates on amplitude levels of the manifold, not on its time axis; the
   consistency of ltep values across geographically separated sites within the same enclosed
   basin confirms this interpretation.

These results motivate future work extending the common manifold to the full global PSMSL
catalogue, investigating the sensitivity of validation skill to the choice of holdout window
(e.g., 1900–1930, 1970–2000), and developing a formal uncertainty quantification framework
for the LTE wavenumber parameter estimates.  Specific high-priority areas for further
investigation include:

- **Honolulu, Hawaii (Site 155)**: Re-optimising the LTE parameters with the wavenumber
  search constrained to the NINO4 basin range (ltep ≈ 10–14) and with the secular trend
  removed from both model and observations, to determine whether Hawaii's sea level genuinely
  tracks the western Pacific warm-pool amplitude mode or whether the Baltic-like ltep ≈ 22
  solution is an artefact of secular confounding in the training record.

- **Fort Denison, Sydney (denison)**: Testing an explicit NINO4-character fit (ltep₁ ≈ 12,
  analogous to the Fremantle solution) to quantify whether the eastern seaboard of Australia
  supports the same Indo-Pacific ENSO amplitude mode that Fremantle carries via the Leeuwin
  Current, and whether this improves or maintains the r_val = 0.411 achieved by the current
  Baltic-template solution.

- **Detrended cross-validation for secular-trend indices** (iode, iodw, tna, tsa): Repeating
  the 1940–1970 holdout analysis after subtracting a linear or quadratic trend from both
  model and observations, to separate the contribution of secular-trend fitting from genuine
  cyclic amplitude-modulation skill.

The mathematical framework underpinning these results is described in *Mathematical Geoenergy*
(Wiley/AGU, 2019) with TeX source at https://github.com/pukpr/gem.

---

*Generated from cross-validation experiments in subdirectories `1` through `tsa`; see each
subdirectory's `lte_results.csv` and `*site1940-1970.png` for individual site detail.  All
correlations quoted are Pearson r computed on the 1940.0–1970.0 decimal-year slice of the
respective `lte_results.csv` (column 2: model output; column 3: observed anomaly).*
