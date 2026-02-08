# GEM-LTE Build and Setup Script
# Automates GNAT installation, build, and deployment

param(
    [switch]$SkipGNAT,
    [switch]$BuildOnly
)

Write-Host "=== GEM-LTE Setup ===" -ForegroundColor Cyan
Write-Host ""

# Check if GNAT is installed
$gnatPath = "C:\GNAT\2021\bin"
$gnatInstalled = Test-Path "$gnatPath\gprbuild.exe"

if (-not $gnatInstalled -and -not $SkipGNAT) {
    Write-Host "GNAT not found at $gnatPath" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Please download and install GNAT Community 2021 from:" -ForegroundColor Yellow
    Write-Host "https://www.adacore.com/download" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "After installation, run this script again." -ForegroundColor Yellow
    
    $response = Read-Host "Open download page in browser? (y/n)"
    if ($response -eq 'y') {
        Start-Process "https://www.adacore.com/download"
    }
    exit 1
}

if ($gnatInstalled) {
    Write-Host "GNAT found at $gnatPath" -ForegroundColor Green
    
    # Add GNAT to PATH
    $env:Path += ";$gnatPath"
    
    # Build the project
    Write-Host ""
    Write-Host "Building lt.exe..." -ForegroundColor Yellow
    gprbuild -P lte.gpr enso_opt
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Build successful" -ForegroundColor Green
        
        # Copy executable to run directory
        Write-Host ""
        Write-Host "Deploying to run directory..." -ForegroundColor Yellow
        
        if (Test-Path "run\lt.exe") {
            Copy-Item "run\lt.exe" "run\lt.exe.bak" -Force
            Copy-Item "experiments\Feb2026\lt.exe" "experiments\Feb2026\lt.exe.bak" -Force
            Write-Host "Backed up previous lt.exe" -ForegroundColor Green
        }
        
        Copy-Item "obj\enso_opt.exe" "run\lt.exe" -Force
        Copy-Item "obj\enso_opt.exe" "experiments\Feb2026\lt.exe" -Force
        Write-Host "Deployed new lt.exe" -ForegroundColor Green
        
        if (-not $BuildOnly) {
            Write-Host ""
            Write-Host "=== Setup Complete ===" -ForegroundColor Cyan
            Write-Host ""
            Write-Host "To run the GUI:" -ForegroundColor Yellow
            Write-Host "  cd experiments\Feb2026" -ForegroundColor White
            Write-Host "  python lte_gui.py" -ForegroundColor White
            Write-Host ""
        }
    } else {
        Write-Host "Build failed" -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "Skipping GNAT check as requested" -ForegroundColor Yellow
}