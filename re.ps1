param (
    [string]$ClimateIndex
)

# Set the environment variable persistently
# [System.Environment]::SetEnvironmentVariable("env:CLIMATE_INDEX", "$ClimateIndex.dat", "User")
$env:CLIMATE_INDEX="$ClimateIndex.dat"

# Delete the file if it exists
$parFile = "lt.exe.$ClimateIndex.par"
if (Test-Path $parFile) {
    Remove-Item $parFile -Force
}

# Copy files
Copy-Item -Path "lt.exe.par.$ClimateIndex" -Destination "lt.exe.par" -Force
Copy-Item -Path "lt.exe.resp.$ClimateIndex" -Destination "lt.exe.resp" -Force