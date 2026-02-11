# GEM-LTE Test Script
# Runs lt.exe and compares last 2 lines of output against reference
# 
# Usage:
#   .\run_test.ps1          # Test with .par files (default)
#   .\run_test.ps1 -json    # Test with .p (JSON) files

param(
    [switch]$json = $false
)

if ($json) {
    Write-Host "Running GEM-LTE test with JSON (.p) files..."
    
    # Ensure JSON files exist
    if (-not (Test-Path "lt.exe.p")) {
        Write-Host "Creating lt.exe.p from lt.exe.par..." -ForegroundColor Yellow
        python par_to_json.py lt.exe.par lt.exe.p
    }
    if (-not (Test-Path "lt.exe.nino4.dat.p")) {
        Write-Host "Creating lt.exe.nino4.dat.p from lt.exe.nino4.dat.par..." -ForegroundColor Yellow
        python par_to_json.py lt.exe.nino4.dat.par lt.exe.nino4.dat.p
    }
    
    # Backup .par files and remove them so JSON is used
    Copy-Item lt.exe.par lt.exe.par.bak -Force
    Copy-Item lt.exe.nino4.dat.par lt.exe.nino4.dat.par.bak -Force
    Remove-Item lt.exe.par -ErrorAction SilentlyContinue
    Remove-Item lt.exe.nino4.dat.par -ErrorAction SilentlyContinue
    
    try {
        # Run lt.exe with JSON files
        .\lt.exe > current_output.txt 2>&1
    } finally {
        # Restore .par files
        Move-Item lt.exe.par.bak lt.exe.par -Force
        Move-Item lt.exe.nino4.dat.par.bak lt.exe.nino4.dat.par -Force
    }
} else {
    Write-Host "Running GEM-LTE test with .par files..."
    
    # Run lt.exe with .par files (standard test)
    .\lt.exe > current_output.txt 2>&1
}

# Get last 2 lines
Get-Content current_output.txt -Tail 2 > current_last2.txt

# Compare with reference
$reference = Get-Content reference_output.txt
$current = Get-Content current_last2.txt

Write-Host "`n=== Reference (last 2 lines) ==="
$reference | ForEach-Object { Write-Host $_ }

Write-Host "`n=== Current Run (last 2 lines) ==="
$current | ForEach-Object { Write-Host $_ }

# Compare
if (Compare-Object $reference $current) {
    Write-Host "`n*** TEST FAILED: Output differs from reference ***" -ForegroundColor Red
    exit 1
} else {
    Write-Host "`n*** TEST PASSED: Output matches reference ***" -ForegroundColor Green
    exit 0
}
