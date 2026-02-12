# Restore parameter files from last passing test commits
# Usage: .\restore_par.ps1

Write-Host 'Restoring parameter files from baseline commits...' -ForegroundColor Cyan

# Restore .par files from commit 18f8709 (test infrastructure baseline)
Write-Host '  Restoring *.par files from 18f8709...' -ForegroundColor Yellow
git checkout 18f8709 -- *.par

# Restore .p (JSON) files from commit f1b989a (JSON test files)  
Write-Host '  Restoring *.p files from f1b989a...' -ForegroundColor Yellow
git checkout f1b989a -- *.p

Write-Host ''
Write-Host 'Done! Parameter files restored.' -ForegroundColor Green
Write-Host 'These should produce cc=0.4107660797 matching reference_output.txt'
Write-Host ''
Write-Host 'Verify with:'
Write-Host '  .\run_test.ps1      # Test with .par files'
Write-Host '  .\run_test.ps1 -j   # Test with JSON .p files'