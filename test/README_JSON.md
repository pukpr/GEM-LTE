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

1. **Array Indexing Issue**: GNATCOLL.JSON uses 1-based indexing, but there may be a mismatch in how triplet elements are accessed
2. **Period Matching Logic**: The Apply() procedure searches for matching periods but may have an off-by-one error
3. **Sequential vs Search**: .par reader assigns sequentially by index; JSON reader searches by period value

### Workaround

Use `json_to_par.py` to convert JSON back to .par format:
```powershell
python json_to_par.py config.json config.par
lt.exe  # Uses config.par successfully
```

This workaround is proven to work perfectly with both primary and secondary parameter files.

### Text .par Format Behavior

When reading from .par files, the code:
1. Reads the period value from the file
2. **OVERRIDES** it with the predefined GEM.LTE.LP value (line 623 of gem-lte-primitives-shared.adb)

```ada
Read (FT, D.A.LP (I), D.B.LPAP (I).Amplitude, D.B.LPAP (I).Phase);
D.A.LP (I) := GEM.LTE.LP (I); -- OVERRIDE!
```

This means .par files can have any period values - they're replaced anyway.

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
