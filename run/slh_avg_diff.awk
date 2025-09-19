awk -F";" '
BEGIN {
    tol = 0.01       # time match tolerance for t - 1
    window = 6       # 13-point window: ±6
}
{
    gsub(/^[ \t]+|[ \t]+$/, "", $2)
    if ($2 == -99999) next
    t = $1 + 0
    v = $2 + 0
    raw_time[n] = t
    raw_value[t] = v
    n++
}
END {
    # Step 1: Perform one-year differencing
    m = 0
    for (i = 0; i < n; i++) {
        t = raw_time[i]
        t_prev = t - 1.0
        found = 0
        for (j = 0; j < n; j++) {
            t_j = raw_time[j]
            if (t_j >= t_prev - tol && t_j <= t_prev + tol) {
                diff_time[m] = t
                diff_value[m] = raw_value[t] - raw_value[t_j]
                m++
                found = 1
                break
            }
        }
        # skip if no match
    }

    if (m < 13) {
        print "Not enough data for 13-point smoothing"
        exit 1
    }

    # Step 2: Apply centered 13-point running mean
    for (i = window; i < m - window; i++) {
        sum = 0
        for (j = i - window; j <= i + window; j++) {
            sum += diff_value[j]
        }
        avg = sum / (2 * window + 1)
        printf "%.4f %.6f\n", diff_time[i], avg
    }
}
' $1
