param (
    [string]$ClimateIdx
)

# Set the environment variable persistently
# [System.Environment]::SetEnvironmentVariable("env:CLIMATE_INDEX", "$ClimateIndex.dat", "User")
$env:CLIMATE_INDEX="$ClimateIdx.dat"

# Delete the file if it exists
$parFile = "lt.exe.$env:CLIMATE_INDEX.par"
if (Test-Path $parFile) {
    notepad $parFile 
}
