# GEM-LTE Test Script
# Runs lt.exe and compares last 2 lines of output against reference
# 
# Usage:
#   .\run_test.ps1          # Test with .par files (default)
#   .\run_test.ps1 -json    # Test with .p (JSON) files

param(
    [switch]$json = $false
)

# Function to check if secondary files are consistent with primary files
function Test-FileConsistency {
    param(
        [string]$extension  # ".par" or ".p"
    )
    
    $primary = "lt.exe$extension"
    $secondary = "lt.exe.nino4.dat$extension"
    
    if (-not (Test-Path $primary)) {
        Write-Host "Warning: Primary file $primary not found" -ForegroundColor Yellow
        return $false
    }
    
    if (-not (Test-Path $secondary)) {
        Write-Host "Warning: Secondary file $secondary not found" -ForegroundColor Yellow
        return $false
    }
    
    $primaryTime = (Get-Item $primary).LastWriteTime
    $secondaryTime = (Get-Item $secondary).LastWriteTime
    
    # Check if files are out of sync (more than 10 seconds apart)
    $timeDiff = [Math]::Abs(($primaryTime - $secondaryTime).TotalSeconds)
    if ($timeDiff -gt 10) {
        Write-Host "Warning: Primary and secondary $extension files have different timestamps" -ForegroundColor Yellow
        Write-Host "  $primary : $primaryTime" -ForegroundColor Yellow
        Write-Host "  $secondary : $secondaryTime" -ForegroundColor Yellow
        Write-Host "  Difference: $([Math]::Round($timeDiff, 1)) seconds" -ForegroundColor Yellow
        
        # Auto-regenerate if needed
        if ($extension -eq ".p") {
            Write-Host "  Auto-regenerating secondary JSON file..." -ForegroundColor Cyan
            python par_to_json.py "lt.exe.nino4.dat.par" $secondary 2>&1 | Out-Null
            if ($LASTEXITCODE -eq 0) {
                Write-Host "  Regenerated successfully" -ForegroundColor Green
            }
        }
        return $false
    }
    
    return $true
}

if ($json) {
    Write-Host "Running GEM-LTE test with JSON (.p) files..."
    
    # Check consistency between .par and .p files
    Write-Host "`nChecking file consistency..."
    Test-FileConsistency ".par" | Out-Null
    Test-FileConsistency ".p" | Out-Null
    
    # Ensure JSON files exist and are up-to-date
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
        # Set test environment variables (TEST_ONLY=TRUE bypasses optimization)
        $env:TEST_ONLY = "TRUE"
        $env:TIMEOUT = "0.0"
        $env:NUMBER_OF_PROCESSORS = "1"
        
        # Run lt.exe with JSON files using -j flag
        .\lt.exe -j > current_output.txt 2>&1
    } finally {
        # Restore .par files
        Move-Item lt.exe.par.bak lt.exe.par -Force
        Move-Item lt.exe.nino4.dat.par.bak lt.exe.nino4.dat.par -Force
    }
} else {
    Write-Host "Running GEM-LTE test with .par files..."
    
    # Check consistency between .par files
    Write-Host "`nChecking file consistency..."
    Test-FileConsistency ".par" | Out-Null
    
    # Set test environment variables (TEST_ONLY=TRUE bypasses optimization)
    $env:TEST_ONLY = "TRUE"
    $env:TIMEOUT = "0.0"
    $env:NUMBER_OF_PROCESSORS = "1"
    
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
