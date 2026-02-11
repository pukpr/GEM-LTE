# Fix Ada style issues (token spacing)
# Adds spaces after commas and around operators

$files = Get-ChildItem src\*.ad? 

foreach ($file in $files) {
    Write-Host "Processing $($file.Name)..."
    $content = Get-Content $file.FullName -Raw
    
    # Fix spacing after commas (but not in strings or comments)
    # Pattern: comma followed by non-whitespace (not inside quotes)
    $lines = Get-Content $file.FullName
    $newLines = @()
    
    foreach ($line in $lines) {
        # Skip comment lines
        if ($line -match '^\s*--') {
            $newLines += $line
            continue
        }
        
        $newLine = $line
        
        # Check if line has comma followed by non-space
        if ($line -match ',\S') {
            # Split on quotes to avoid modifying strings
            $parts = $line -split '(".*?")', -1, 'RegexMatch'
            for ($i = 0; $i -lt $parts.Length; $i += 2) {
                # Only modify non-quoted parts (even indices)
                $parts[$i] = $parts[$i] -replace ',([^\s\)])', ', $1'
            }
            $newLine = $parts -join ''
        }
        
        $newLines += $newLine
    }
    
    # Write back
    $newLines | Set-Content $file.FullName
}

Write-Host "`nDone! Rebuild to check for remaining issues."
