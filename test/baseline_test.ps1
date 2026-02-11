# Baseline Test - Capture iteration trace with fixed seed
# This establishes the baseline behavior before refactoring

param(
    [int]$Seed = 42,
    [int]$MaxLoops = 500,
    [string]$OutputFile = "baseline_trace.txt"
)

$env:SEED = "$Seed"
$env:MAXLOOPS = "$MaxLoops"

Write-Host "Running baseline test with:"
Write-Host "  SEED = $Seed"
Write-Host "  MAXLOOPS = $MaxLoops"
Write-Host "  Output = $OutputFile"
Write-Host ""

.\lt.exe > $OutputFile 2>&1

# Extract key metrics
$ccLine = Get-Content $OutputFile | Select-String "^cc " | Select-Object -Last 1
$totalLines = (Get-Content $OutputFile).Count

Write-Host "Baseline Results:"
Write-Host "  Total output lines: $totalLines"
Write-Host "  Final: $($ccLine.Line)"
Write-Host "  Output saved to: $OutputFile"
Write-Host ""

# Create summary
$summary = @"
Baseline Test Summary
=====================
Date: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")
Seed: $Seed
MaxLoops: $MaxLoops
Total Lines: $totalLines
Final Result: $($ccLine.Line)

This file captures the exact behavior before refactoring.
Use this to verify refactored code produces identical results.

To compare after refactoring:
  .\baseline_test.ps1 -Seed $Seed -MaxLoops $MaxLoops -OutputFile refactored_trace.txt
  Compare-Object (Get-Content baseline_trace.txt) (Get-Content refactored_trace.txt)
"@

$summary | Out-File "baseline_summary.txt"
Write-Host "Summary saved to: baseline_summary.txt"
