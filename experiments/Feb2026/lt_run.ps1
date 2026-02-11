param (
    [string]$Index,
    [string]$Time
)

Write-Host "Processing file: $Index"

$env:METRIC="cc"
$env:TIMEOUT="$Time"
$env:TRAIN_START="1940"
$env:TRAIN_STOP="1970"
$env:CLIMATE_INDEX="$Index.dat"
$env:IDATE="1920.9"

del lt.exe.$Index.dat.par

..\lt

python3  ..\plot.py $Index Feb2026 $env:TRAIN_START $env:TRAIN_STOP  0




