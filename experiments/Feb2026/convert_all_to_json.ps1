#!/usr/bin/env pwsh
# Convert all *.par files in Feb2026 site directories to JSON format

$ErrorActionPreference = "Stop"
$script_dir = Split-Path -Parent $MyInvocation.MyCommand.Path
$converter = Join-Path $script_dir "..\..\test\par_to_json.py"

if (-not (Test-Path $converter)) {
    Write-Error "Cannot find converter: $converter"
    exit 1
}

Write-Host "Converting all *.par files to JSON format..."
Write-Host "Converter: $converter"
Write-Host ""

# Get all subdirectories (site directories)
$dirs = Get-ChildItem -Path $script_dir -Directory | 
        Where-Object { $_.Name -notin @("locs", "scripts", "rlr_data") }

$total_converted = 0
$failed = @()

foreach ($dir in $dirs) {
    $par_files = Get-ChildItem -Path $dir.FullName -Filter "*.par" -File
    
    if ($par_files.Count -eq 0) {
        continue
    }
    
    Write-Host "Processing directory: $($dir.Name)"
    
    foreach ($par in $par_files) {
        $json_file = $par.FullName -replace '\.par$', '.p'
        
        Write-Host "  $($par.Name) -> $([System.IO.Path]::GetFileName($json_file))"
        
        try {
            $result = & python $converter $par.FullName $json_file 2>&1
            if ($LASTEXITCODE -ne 0) {
                Write-Warning "    FAILED: $result"
                $failed += "$($dir.Name)/$($par.Name)"
            } else {
                $total_converted++
            }
        } catch {
            Write-Warning "    ERROR: $_"
            $failed += "$($dir.Name)/$($par.Name)"
        }
    }
}

Write-Host ""
Write-Host "========================================="
Write-Host "Conversion complete!"
Write-Host "  Converted: $total_converted files"
Write-Host "  Failed:    $($failed.Count) files"
Write-Host "========================================="

if ($failed.Count -gt 0) {
    Write-Host ""
    Write-Host "Failed files:"
    foreach ($f in $failed) {
        Write-Host "  $f"
    }
    exit 1
}

exit 0
