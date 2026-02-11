# GEM-LTE Test Script
# Runs lt.exe and compares last 2 lines of output against reference

Write-Host "Running GEM-LTE test..."

# Run lt.exe and capture output
.\lt.exe > current_output.txt 2>&1

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
