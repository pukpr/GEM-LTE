# GEM-LTE Test Directory

This directory contains a simple regression test for the GEM-LTE climate modeling code.

## Test Setup

The test uses the `enso_opt` executable (copied as `lt.exe`) with a single-pass configuration to verify consistent model behavior.

### Files

- `lt.exe` - Copy of `enso_opt.exe` executable
- `lt.exe.par` - Model parameters
- `lt.exe.resp` - Configuration settings (includes TIMEOUT=0.0, TRIGGER=0.0)
- `*.dat` - Input data files (nino4.dat, dlod3.dat, etc.)
- `run_test.ps1` - Test script
- `reference_output.txt` - Reference output (last 2 lines)

## Running the Test

```powershell
.\run_test.ps1
```

The script:
1. Runs `lt.exe` with the configured parameters
2. Captures the last 2 lines of output
3. Compares against the reference output
4. Returns exit code 0 (pass) or 1 (fail)

## Creating Reference Output

To update the reference output after code changes:

```powershell
.\lt.exe > test_output.txt 2>&1
Get-Content test_output.txt -Tail 2 > reference_output.txt
```

## Notes

- The test uses `nino4.dat` as the climate index (configured in lt.exe.resp)
- Configuration files set TIMEOUT=0.0 and TRIGGER=0.0 for deterministic single-pass execution
- The last two lines show the final dLOD correlation and year values
- Minor numerical differences may occur due to optimization iteration variations
