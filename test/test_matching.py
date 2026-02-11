#!/usr/bin/env python3
"""Test period matching logic"""

# Predefined LP periods from GEM.LTE.LP array
LP_PERIODS = [
    27.32166155400,
    27.21222081500,
    1095.17502233364,
    9.08286783775,
    6.81670784238,
    13.63341568476,
    13.66083077700,
    13.60611040750,
    27.55454988600,
    13.77727494300,
    6793.47648813065,
    1616.30271425126,
    27.09267692660,
    9.09502787633,
    13.69115772864,
    1305.66025845931,
    27.44323926226,
    27.66676714572,
    2190.35004466729,
    26.98505934729,
    6167.20701373239,
    -3232.60542850251,
    -2120.90900683440,
    27.33282433123,
    13.71880577217,
    9.13295078376,
    9.12068919638,
    9.10846048884,
    3396.73824406533,
]

# Periods from .par file (lines 19-47)
PAR_PERIODS = [
    27.32166155400,
    27.21222081500,
    1095.17502233364,
    9.08286783775,
    6.81670784238,
    13.63341568476,
    13.66083077700,
    13.60611040750,
    27.55454988600,
    13.77727494300,
    6793.47648813065,
    1616.30271425126,
    27.09267692660,
    9.09502787633,
    13.69115772864,
    1305.66025845931,
    27.44323926226,
    27.66676714572,
    2190.35004466729,
    26.98505934729,
    6167.20701373239,
    -3232.60542850251,
    -2120.90900683440,
    27.33282433123,
    13.71880577217,
    9.13295078376,
    9.12068919638,
    9.10846048884,
    3396.73824406533,
]

def test_matching():
    """Test period matching with 1% tolerance"""
    tolerance = 0.01
    
    print("Testing period matching:")
    print("=" * 60)
    
    for i, par_period in enumerate(PAR_PERIODS):
        # Find best match
        best_match = None
        best_idx = None
        best_error = float('inf')
        
        for j, lp_period in enumerate(LP_PERIODS):
            if abs(lp_period) > 0.0:
                error = abs(par_period - lp_period)
                threshold = abs(lp_period) * tolerance
                if error <= threshold and error < best_error:
                    best_match = lp_period
                    best_idx = j
                    best_error = error
        
        if best_match is not None:
            if best_idx != i:
                print(f"{i+1:2}. PAR {par_period:20.11f} -> LP[{best_idx+1}] "
                      f"{best_match:20.11f} (error {best_error:.2e})")
            # else: exact match at same index, no output needed
        else:
            print(f"{i+1:2}. PAR {par_period:20.11f} -> NO MATCH!")

if __name__ == "__main__":
    test_matching()
