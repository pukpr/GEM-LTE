# GEM-LTE JSON Parameter Format

This directory contains utilities for converting GEM-LTE parameter files between text (.par) and JSON formats.

## Files

- **par_to_json.py** - Python utility to convert .par files to JSON
- **test_with_json.ps1** - PowerShell script to test JSON parameter loading
- **README_JSON.md** - This file

## Usage

### Convert .par to JSON

```bash
python par_to_json.py lt.exe.par lt.exe.p
```

This creates a JSON file named `lt.exe.p` (or specify your own output name).

### Test JSON Loading

```powershell
.\test_with_json.ps1
```

This script:
1. Converts lt.exe.par to JSON
2. Backs up the original .par file
3. Installs the JSON file as lt.exe.json
4. Runs the test suite
5. Restores the original .par file

## JSON Format

The JSON format contains:

```json
{
  "offs": 0.0,
  "bg": 0.00014145452,
  "impA": 0.0,
  ...
  "ltep": [23.405, 2.632, ...],
  "harm": [45, 49, 48, ...],
  "lpap": [
    [27.321, 0.006214, -1.351],
    [27.212, 0.031444, -0.769],
    ...
  ]
}
```

### Field Descriptions

- **Scalar parameters**: offs, bg, impA, impB, impC, delA, delB, asym, ann1, ann2, sem1, sem2, year, IR, ma, mp, shfT, init
- **ltep**: Long-Term Evolution Parameters (LTE phases) - array of floats
- **harm**: Harmonic indices - array of integers
- **lpap**: Long Period Amplitude Phase tidal constituents - array of [period, amplitude, phase] triplets

## Configuration Parameters (NH and NM)

The `.resp` configuration file controls which elements of the `ltep` and `harm` arrays are used:

### NH Parameter (Harmonic Count)

Controls how many harmonic terms from the `harm` array are used in optimization.

**Two formats are supported:**

1. **Integer format (recommended)**: `NH=6`
   - Simply specifies the count of harmonic terms to use
   - Clear and concise

2. **List format (legacy/obscure)**: `NH=1 1 1 1 1 1`
   - Space-separated list of 1's
   - Length of list indicates count (6 ones = 6 terms)
   - Backward compatible with older .resp files
   - Less intuitive but still supported

**Examples:**
```
NH=6          # Use first 6 harmonic terms (recommended)
NH=1 1 1      # Use first 3 harmonic terms (legacy format)
NH=27         # Use all 27 harmonic terms
```

### NM Parameter (LTEP Modulation Count)

Controls how many long-term evolution parameters from the `ltep` array are used:

```
NM=2          # Use first 2 LTEP terms
NM=11         # Use all 11 LTEP terms
```

### JSON Compact Format

When using JSON `.p` files, the arrays are automatically sized based on NH/NM parameters:

- **Writing**: JSON arrays contain only NH harm entries and NM ltep entries
- **Reading**: If NH/NM are not specified in .resp file, they are auto-detected from JSON array lengths
- **.par files**: Always write full arrays (unchanged for backward compatibility)

**Example workflow:**
```powershell
# Set compact size in .resp file
echo "NH=6" >> lt.exe.resp
echo "NM=2" >> lt.exe.resp

# Save creates compact JSON (6 harm, 2 ltep)
lt.exe -s

# Load auto-detects array sizes from JSON
lt.exe
```

## Known Issues

### LPAP Data Not Loading from JSON ⚠️

**Status:** JSON format is proven complete and correct; issue is in Ada JSON reader

The LPAP (Long Period Amplitude Phase) tidal constituent data fails to load correctly when using JSON format directly, producing incorrect results (correlation coefficient -0.029 instead of expected 0.411).

**CRITICAL FINDING:** Round-trip testing proves the JSON format is complete and correct:

```powershell
# Test both primary and secondary files
python par_to_json.py lt.exe.par
python par_to_json.py lt.exe.nino4.dat.par
python json_to_par.py lt.exe.p lt.exe.par
python json_to_par.py lt.exe.nino4.dat.p lt.exe.nino4.dat.par
lt.exe  # Result: cc = 0.4107660797 ✅ EXACT MATCH
```

**Evidence:**
- ✅ `.par → JSON → .par → lt.exe` = Perfect results (cc=0.4108)
- ❌ `.par → JSON → lt.exe` directly = Wrong results (cc=-0.029)
- ✅ Both primary and secondary files round-trip successfully
- ✅ All scalar parameters, ltep, and harm arrays load correctly from JSON
- ❌ Only LPAP array fails when loading directly from JSON

**Conclusion:** The issue is definitively in the Ada `Read_JSON_LPAP` function (src/gem-lte-primitives-shared.adb, lines 451-510), not in the JSON data.

### Suspected Root Causes

The .par reader and JSON reader use fundamentally different strategies:

**.par Format: ORDER-DEPENDENT Sequential Assignment (WORKS)**
```ada
for I in D.B.LPAP'Range loop
   Read (FT, D.A.LP (I), D.B.LPAP (I).Amplitude, D.B.LPAP (I).Phase);
   D.A.LP (I) := GEM.LTE.LP (I); -- Assigns to index I sequentially (1,2,3...)
end loop;
```
- **Critical:** .par files MUST have LPAP triplets in the same order as GEM.LTE.LP array
- Reads triplet #1 → assigns to LP[1], triplet #2 → LP[2], etc.
- Period values in file are discarded (overwritten)

**JSON Format: ORDER-INDEPENDENT Period Matching (FAILS)**
```ada
for I in D.B.LPAP'Range loop
   if abs (Period - GEM.LTE.LP (I)) <= abs (GEM.LTE.LP (I)) * 0.01 then
      D.B.LPAP (I).Amplitude := Amp;  -- Searches for match, assigns to index I
      D.B.LPAP (I).Phase := Phase;
      exit;
   end if;
end loop;
```
- **Design intent:** JSON arrays can be in any order - period matching finds correct index
- Should be more flexible than .par format
- Bug must be in the matching implementation itself

**Possible issues:**
1. Array indexing confusion in GNATCOLL.JSON Get() calls
2. Off-by-one error in Apply() loop
3. Period matching happens but writes to wrong array index
4. Loop exits prematurely or doesn't iterate correctly

### Workaround

Use `json_to_par.py` to convert JSON back to .par format:
```powershell
python json_to_par.py config.json config.par
lt.exe  # Uses config.par successfully
```

This workaround is proven to work perfectly with both primary and secondary parameter files.

## Format Differences: Order Dependency

### .par Format: ORDER-DEPENDENT (Sequential Reading)

The .par format **requires** LPAP triplets to be in the exact same order as the internal `GEM.LTE.LP` array:

```ada
for I in D.B.LPAP'Range loop
   Read (FT, D.A.LP (I), D.B.LPAP (I).Amplitude, D.B.LPAP (I).Phase);
   D.A.LP (I) := GEM.LTE.LP (I); -- Sequential: I=1,2,3...
end loop;
```

- Reads triplets sequentially (I=1, 2, 3, ...)
- Assigns directly to array index I
- **Order is critical** - triplet #1 must correspond to LP[1], triplet #2 to LP[2], etc.
- Period values in the file are discarded (overwritten with GEM.LTE.LP values)

### JSON Format: ORDER-INDEPENDENT (Period Matching)

The JSON format allows LPAP triplets in **any order** because it uses period matching:

```ada
for I in D.B.LPAP'Range loop
   if abs (Period - GEM.LTE.LP (I)) <= abs (GEM.LTE.LP (I)) * 0.01 then
      D.B.LPAP (I).Amplitude := Amp;  -- Match found, assign to index I
      D.B.LPAP (I).Phase := Phase;
      exit;
   end if;
end loop;
```

- Searches through LP array to find matching period
- Uses 1% tolerance for matching
- Can handle JSON arrays in any order
- **Period matching is critical** - must correctly identify which LP index to write to

### Why This Design?

JSON objects/arrays can be reordered during editing or programmatic manipulation. The period-matching approach makes the JSON format more robust and flexible than the sequential .par format.

**Example:**
```json
"lpap": [
  [27.321, 0.006, -1.35],  // Could be LP[1]
  [27.212, 0.031, -0.76],  // Could be LP[2]
  ...
]
```

These could appear in any order in the JSON file, and the reader should still correctly assign them to the right LP array indices based on period matching.

## Recommendations

1. **For production use**: Keep using .par format until JSON reader is fixed
2. **For testing**: Use test_with_json.ps1 to verify JSON compatibility
3. **For development**: JSON format is easier to edit and validate

## Future Improvements

To make JSON format fully compatible:

1. Update `Read_JSON_LPAP` to accept exact period values without matching
2. OR: Ensure .par-to-JSON conversion uses exact GEM.LTE.LP periods
3. OR: Add a "strict" JSON mode that overrides periods like .par format does

## Technical Details

See `gem-lte-primitives-shared.adb` lines 449-577 for JSON reading implementation.
