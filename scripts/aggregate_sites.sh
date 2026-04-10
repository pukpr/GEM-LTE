#!/bin/bash

OUTPUT="filtered_psmsl_1880_2026.txt"
> "$OUTPUT"   # empty the output file

for N in $(seq 1 999); do
    URL="https://psmsl.org/data/obtaining/rlr.monthly.data/${N}.rlrdata"

    echo "Fetching $URL"

    # Download file silently
    DATA=$(curl -s "$URL")

    # Skip if file not found or empty
    if [[ -z "$DATA" ]]; then
        continue
    fi

    # Extract rows where first column is between 2020.0 and 2025.0
    echo "$DATA" | awk -F';' -v N="$N"  '
       /^[[:space:]]*[0-9]/ {
           gsub(/^[[:space:]]+/, "", $1)   # remove leading spaces from field 1
           date = $1 + 0                   # convert to float
           if (date >= 1880.0 && date <= 2026.0)
               print N " " date " " $2
       }
    ' >> "$OUTPUT"
done

echo "Done. Results saved to $OUTPUT"

