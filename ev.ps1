# Get all environment variables
$allVariables = Get-ChildItem Env:

# Get system environment variables
$systemVariables = [System.Environment]::GetEnvironmentVariables("Machine")

# Define a list of common environment variable names to exclude
$commonVarsToExclude = @("ProgramFiles", "ProgramData", "SystemRoot", "windir", "COMPUTERNAME", "USERDOMAIN", 
"ALLUSERSPROFILE",
"CommonProgramFiles",
"CommonProgramFiles(x86)",
"CommonProgramW6432",
"LOGONSERVER",
"OMSYSHOME",
"OneDrive",
"OneDriveConsumer",
"ProgramFiles(x86)",
"ProgramW6432",
"SESSIONNAME",
"SystemDrive",
"USERDOMAIN_ROAMINGPROFILE",
"fm7home",
"USERNAME", "USERPROFILE", "HOMEPATH", "HOMEDRIVE", "APPDATA", "LOCALAPPDATA", "PUBLIC", "TEMP", "TMP")

# Find user-defined variables by excluding system variables and common variables
$userVariables = $allVariables | Where-Object { -not $systemVariables.ContainsKey($_.Name) -and $commonVarsToExclude -notcontains $_.Name }

# Display user-defined variables
$userVariables

