# JSON Parameter Support - Experiments Directory

## Overview
All experiment parameter files in `experiments/Feb2026/` have been converted to JSON format and the GUI launcher has been updated with a JSON checkbox control.

## What Was Done

### 1. Bulk Conversion (159 Files)
Created `convert_all_to_json.ps1` utility that:
- Scans all site directories in Feb2026/
- Converts all `*.par` files to `*.p` JSON files
- Successfully converted **159 files** across **83 site directories**
- Zero failures

**Converted files:**
- Primary: `lt.exe.p` (includes LPAP tidal constituents)
- Secondary: `lt.exe.<index>.dat.p` (index-specific parameters)

### 2. GUI Update - JSON Checkbox
Modified `lte_gui.py` to add JSON loading control:

**UI Changes:**
- Added **"JSON" checkbox** on the "Run lt" line
- Located immediately after "Run lt" button
- **Defaults to CHECKED** (JSON mode enabled)

**Code Changes:**
- `resolve_lt_cmd()` now accepts `use_json` parameter
- Checkbox state controls whether `-j` flag is passed to `lt.exe`
- When checked: runs `lt.exe -j` (JSON-only mode)
- When unchecked: runs `lt.exe` (standard .par mode)

## Usage

### Running Experiments from GUI
1. Launch `python lte_gui.py` from Feb2026 directory
2. Select a site directory from the list
3. Configure timeout and test interval
4. **JSON checkbox:**
   - ✅ **CHECKED** (default): Load from `.p` JSON files
   - ☐ **UNCHECKED**: Load from `.par` text files
5. Click "Run lt"

### Re-converting Files
If parameters are updated and you need to regenerate JSON files:

```powershell
cd experiments/Feb2026
.\convert_all_to_json.ps1
```

This will update all JSON files to match the current `.par` files.

## Benefits

1. **Programmatic Access**: JSON files are easier to parse and modify with scripts
2. **Version Control**: JSON diffs are more readable than .par format
3. **Validation**: JSON structure can be schema-validated
4. **Flexibility**: Users can choose format per run via checkbox
5. **Auto-Sync**: Optimization updates both .par and .p files automatically
6. **Backward Compatible**: Unchecking JSON box uses legacy .par files

## Files Modified

### Created:
- `experiments/Feb2026/convert_all_to_json.ps1` - Bulk conversion utility
- `experiments/Feb2026/*/lt.exe.p` - Primary JSON files (83 files)
- `experiments/Feb2026/*/lt.exe.*.dat.p` - Secondary JSON files (76 files)

### Modified:
- `experiments/Feb2026/lte_gui.py` - Added JSON checkbox control
- `experiments/Feb2026/lt.exe` - Updated binary with latest JSON support

## Architecture

### JSON File Structure
```json
{
  "offs": 0.0,
  "bg": 0.00014145452,
  "impA": 0.0,
  ...
  "lpap": [
    [27.321582, 0.00123, 1.234],  // [period, amplitude, phase]
    [13.660791, 0.00456, 2.345],
    ...
  ],
  "ltep": [...],
  "harm": [1, 2, 3, 4, 0]  // terminated by 0
}
```

### Dual Format Support
Both `.par` and `.p` files are maintained:
- `.par` files: Legacy text format, human-readable
- `.p` files: JSON format, machine-readable
- Kept synchronized by `Save()` procedure
- GUI checkbox selects which format to load

## Testing
All 159 conversions completed successfully with zero failures.
JSON files have been validated against round-trip testing in the test suite.

## Commits
- `240d0a5`: Add JSON file writing to Save procedure
- `cca8a15`: Convert all experiment *.par files to JSON and add GUI checkbox
