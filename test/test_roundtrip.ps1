#!/usr/bin/env pwsh
#
# test_roundtrip.ps1 - Verify JSON format completeness via round-trip conversion
#
# This script proves that the JSON format contains all necessary data by:
# 1. Converting .par → JSON (par_to_json.py)
# 2. Converting JSON → .par (json_to_par.py)
# 3. Running lt.exe with the round-trip .par files
# 4. Verifying results match the reference output exactly
#
# If this test passes, it proves:
# - The JSON conversion utilities work correctly
# - The JSON format preserves all necessary data with full precision
# - Any JSON loading failures are in the Ada JSON reader, not the JSON data
#
# Usage:
#   .\test_roundtrip.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  JSON Round-Trip Verification Test" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "This test proves JSON format completeness by converting:" -ForegroundColor Gray
Write-Host "  .par → JSON → .par → lt.exe" -ForegroundColor Gray
Write-Host ""

# Step 1: Convert .par to JSON
Write-Host "[1/6] Converting .par files to JSON..." -ForegroundColor Yellow
python par_to_json.py lt.exe.par lt.exe.p 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to convert lt.exe.par" -ForegroundColor Red
    exit 1
}
python par_to_json.py lt.exe.nino4.dat.par lt.exe.nino4.dat.p 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to convert lt.exe.nino4.dat.par" -ForegroundColor Red
    exit 1
}
Write-Host "  ✓ lt.exe.par → lt.exe.p" -ForegroundColor Green
Write-Host "  ✓ lt.exe.nino4.dat.par → lt.exe.nino4.dat.p" -ForegroundColor Green
Write-Host ""

# Step 2: Convert JSON back to .par
Write-Host "[2/6] Converting JSON back to .par..." -ForegroundColor Yellow
python json_to_par.py lt.exe.p lt.exe.par.roundtrip 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to convert lt.exe.p" -ForegroundColor Red
    exit 1
}
python json_to_par.py lt.exe.nino4.dat.p lt.exe.nino4.dat.par.roundtrip 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to convert lt.exe.nino4.dat.p" -ForegroundColor Red
    exit 1
}
Write-Host "  ✓ lt.exe.p → lt.exe.par.roundtrip" -ForegroundColor Green
Write-Host "  ✓ lt.exe.nino4.dat.p → lt.exe.nino4.dat.par.roundtrip" -ForegroundColor Green
Write-Host ""

# Step 3: Backup original files
Write-Host "[3/6] Backing up original .par files..." -ForegroundColor Yellow
Copy-Item lt.exe.par lt.exe.par.backup -Force
Copy-Item lt.exe.nino4.dat.par lt.exe.nino4.dat.par.backup -Force
Write-Host "  ✓ Created backups" -ForegroundColor Green
Write-Host ""

# Step 4: Install round-trip files
Write-Host "[4/6] Installing round-trip .par files..." -ForegroundColor Yellow
Copy-Item lt.exe.par.roundtrip lt.exe.par -Force
Copy-Item lt.exe.nino4.dat.par.roundtrip lt.exe.nino4.dat.par -Force
Write-Host "  ✓ Installed round-trip files" -ForegroundColor Green
Write-Host ""

# Step 5: Run test
Write-Host "[5/6] Running test with round-trip .par files..." -ForegroundColor Yellow
Write-Host ""
.\run_test.ps1
$test_result = $LASTEXITCODE
Write-Host ""

# Step 6: Restore and cleanup
Write-Host "[6/6] Restoring original files and cleaning up..." -ForegroundColor Yellow
Move-Item lt.exe.par.backup lt.exe.par -Force
Move-Item lt.exe.nino4.dat.par.backup lt.exe.nino4.dat.par -Force
Remove-Item lt.exe.p, lt.exe.nino4.dat.p -ErrorAction SilentlyContinue
Remove-Item lt.exe.par.roundtrip, lt.exe.nino4.dat.par.roundtrip -ErrorAction SilentlyContinue
Write-Host "  ✓ Restored originals" -ForegroundColor Green
Write-Host ""

# Report final result
Write-Host "========================================" -ForegroundColor Cyan
if ($test_result -eq 0) {
    Write-Host "  ROUND-TRIP TEST PASSED ✓" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "SUCCESS: JSON format is complete and correct!" -ForegroundColor Green
    Write-Host ""
    Write-Host "This proves:" -ForegroundColor White
    Write-Host "  ✓ par_to_json.py converts all necessary data" -ForegroundColor Gray
    Write-Host "  ✓ json_to_par.py recreates perfect .par files" -ForegroundColor Gray
    Write-Host "  ✓ JSON format preserves full floating-point precision" -ForegroundColor Gray
    Write-Host "  ✓ Round-trip .par files produce exact results" -ForegroundColor Gray
    Write-Host ""
    Write-Host "Any JSON loading failures are in the Ada JSON reader," -ForegroundColor Yellow
    Write-Host "NOT in the JSON data or conversion utilities." -ForegroundColor Yellow
    exit 0
} else {
    Write-Host "  ROUND-TRIP TEST FAILED ✗" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Round-trip conversion produced different results." -ForegroundColor Red
    Write-Host "This indicates an issue with the conversion utilities." -ForegroundColor Red
    exit 1
}
