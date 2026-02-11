# Master setup script for GEM-LTE
# Detects available options and guides user

Write-Host @"
===========================================
   GEM-LTE Setup
   GeoEnergyMath Laplace's Tidal Equation
===========================================
"@ -ForegroundColor Cyan

Write-Host ""

# Check if pre-built executable exists
$hasExe = Test-Path "run\lt.exe"

if ($hasExe) {
    Write-Host "✓ Pre-built lt.exe found" -ForegroundColor Green
    Write-Host ""
    Write-Host "You can run GEM-LTE immediately!" -ForegroundColor Green
    Write-Host ""
    Write-Host "Options:" -ForegroundColor Yellow
    Write-Host "  [1] Launch GUI (requires Python)" -ForegroundColor White
    Write-Host "  [2] Rebuild from source (requires GNAT)" -ForegroundColor White
    Write-Host "  [3] Exit" -ForegroundColor White
    Write-Host ""
    
    $choice = Read-Host "Select option (1-3)"
    
    switch ($choice) {
        "1" {
            # Check Python
            if (Get-Command python -ErrorAction SilentlyContinue) {
                Write-Host "Launching GUI..." -ForegroundColor Green
                Set-Location experiments\Feb2026
                python lte_gui.py
            } else {
                Write-Host "ERROR: Python not found" -ForegroundColor Red
                Write-Host "Please install Python from https://www.python.org" -ForegroundColor Yellow
            }
        }
        "2" {
            .\setup_with_gnat.ps1
        }
        "3" {
            exit 0
        }
        default {
            Write-Host "Invalid option" -ForegroundColor Red
        }
    }
} else {
    Write-Host "Pre-built executable not found." -ForegroundColor Yellow
    Write-Host "You need to build from source." -ForegroundColor Yellow
    Write-Host ""
    .\setup_with_gnat.ps1
}