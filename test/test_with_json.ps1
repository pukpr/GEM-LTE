#!/usr/bin/env pwsh
#
# test_with_json.ps1 - Convert .par to JSON and run test
#
# This script automates the process of:
# 1. Converting lt.exe.par to lt.exe.p (JSON format)
# 2. Renaming files to use JSON format
# 3. Running the test suite
# 4. Restoring original files
#
# NOTE: JSON format currently has limitations:
#   - LPAP (tidal constituent) matching uses period values to find
#     the correct entry in the predefined GEM.LTE.LP array
#   - Periods must match within 1% tolerance
#   - If periods don't match, LPAP data is not loaded properly
#   - This can cause test failures due to missing tidal constituents
#
# For production use, consider keeping .par format or fixing the
# JSON reader to handle exact period values.
#
# Usage:
#   .\test_with_json.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  GEM-LTE JSON Parameter Test" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Step 1: Convert .par files to JSON
Write-Host "[1/6] Converting .par files to JSON..." -ForegroundColor Yellow
python par_to_json.py lt.exe.par lt.exe.p
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Conversion of lt.exe.par failed" -ForegroundColor Red
    exit 1
}
python par_to_json.py lt.exe.nino4.dat.par lt.exe.nino4.dat.p
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Conversion of lt.exe.nino4.dat.par failed" -ForegroundColor Red
    exit 1
}
Write-Host ""

# Step 2: Backup original files
Write-Host "[2/6] Backing up original files..." -ForegroundColor Yellow
if (Test-Path "lt.exe.par.backup") {
    Remove-Item "lt.exe.par.backup"
}
if (Test-Path "lt.exe.nino4.dat.par.backup") {
    Remove-Item "lt.exe.nino4.dat.par.backup"
}
Copy-Item "lt.exe.par" "lt.exe.par.backup"
Copy-Item "lt.exe.nino4.dat.par" "lt.exe.nino4.dat.par.backup"
Write-Host "  - Backed up lt.exe.par -> lt.exe.par.backup" -ForegroundColor Gray
Write-Host "  - Backed up lt.exe.nino4.dat.par -> lt.exe.nino4.dat.par.backup" -ForegroundColor Gray
Write-Host ""

# Step 3: Move JSON files into place
Write-Host "[3/6] Installing JSON parameter files..." -ForegroundColor Yellow
if (Test-Path "lt.exe.json") {
    Remove-Item "lt.exe.json"
}
if (Test-Path "lt.exe.nino4.dat.json") {
    Remove-Item "lt.exe.nino4.dat.json"
}
Move-Item "lt.exe.p" "lt.exe.json"
Move-Item "lt.exe.nino4.dat.p" "lt.exe.nino4.dat.json"
Write-Host "  - Installed lt.exe.p -> lt.exe.json" -ForegroundColor Gray
Write-Host "  - Installed lt.exe.nino4.dat.p -> lt.exe.nino4.dat.json" -ForegroundColor Gray
Write-Host ""

# Step 4: Run test
Write-Host "[4/6] Running test with JSON parameters..." -ForegroundColor Yellow
Write-Host ""
.\run_test.ps1
$test_result = $LASTEXITCODE
Write-Host ""

# Step 5: Restore original files
Write-Host "[5/6] Restoring original files..." -ForegroundColor Yellow
if (Test-Path "lt.exe.json") {
    Remove-Item "lt.exe.json"
}
if (Test-Path "lt.exe.nino4.dat.json") {
    Remove-Item "lt.exe.nino4.dat.json"
}
Move-Item "lt.exe.par.backup" "lt.exe.par" -Force
Move-Item "lt.exe.nino4.dat.par.backup" "lt.exe.nino4.dat.par" -Force
Write-Host "  - Restored lt.exe.par" -ForegroundColor Gray
Write-Host "  - Restored lt.exe.nino4.dat.par" -ForegroundColor Gray
Write-Host ""

# Report final result
Write-Host "========================================" -ForegroundColor Cyan
if ($test_result -eq 0) {
    Write-Host "  JSON TEST PASSED ✓" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "JSON parameter format is working correctly!" -ForegroundColor Green
    Write-Host "The lt.exe.par file can be converted to JSON format." -ForegroundColor Green
    exit 0
} else {
    Write-Host "  JSON TEST FAILED ✗" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "JSON parameter loading produced different results." -ForegroundColor Red
    Write-Host "This may indicate a conversion issue or JSON reader bug." -ForegroundColor Red
    exit 1
}
