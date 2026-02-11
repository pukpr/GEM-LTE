#!/usr/bin/env python3
"""
par_to_json.py - Convert GEM-LTE .par files to JSON format

Converts text-based .par parameter files to JSON format for easier
configuration management and programmatic access.

IMPORTANT: Order and Period Matching
=====================================
.par format: LPAP triplets MUST be in sequential order (1,2,3...)
             matching the internal GEM.LTE.LP array indices.

JSON format: LPAP triplets CAN be in any order because the Ada JSON
             reader uses period matching to find the correct LP index.

This converter preserves the sequential order from .par files, but the
predefined LP periods are used to ensure exact matching with the Ada
JSON reader's period-matching logic (1% tolerance).

Usage:
    python par_to_json.py <input.par> [output.json]
    
If output filename is not specified, uses input basename with .p extension.

Format:
    .par files contain tagged values:
        offs    0.12345
        impA    1.23456
        ltep    2.71828  3.14159  1.41421
        harm    45
    
    JSON output:
        {
            "offs": 0.12345,
            "impA": 1.23456,
            "ltep": [2.71828, 3.14159, 1.41421],
            "harm": [45, ...]
        }
"""

import sys
import json
import re
from pathlib import Path
from typing import Dict, List, Union, Optional


# Predefined LP periods from GEM.LTE.LP array (computed from Doodson arguments)
# These are the canonical periods used by the Ada code
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


def match_period(file_period: float, tolerance: float = 0.01) -> Optional[float]:
    """
    Match a period from the .par file to a predefined LP period.
    
    Since .par file periods are identical to LP periods, this function
    looks for an exact match first, then falls back to tolerance matching.
    
    Args:
        file_period: Period value from .par file
        tolerance: Matching tolerance (default 1% = 0.01)
        
    Returns:
        Matched predefined period, or None if no match found
    """
    # First try exact match (periods in .par are identical to LP)
    for lp_period in LP_PERIODS:
        if file_period == lp_period:
            return lp_period
    
    # Fall back to tolerance matching for robustness
    for lp_period in LP_PERIODS:
        if abs(lp_period) > 0.0:
            if abs(file_period - lp_period) <= abs(lp_period) * tolerance:
                return lp_period
    
    return None


def parse_par_file(filename: str) -> Dict[str, Union[float, List[float], List[int]]]:
    """
    Parse a .par file and return a dictionary suitable for JSON serialization.
    
    Uses predefined LP periods for LPAP triplets to ensure exact matching
    with Ada JSON reader's period-matching logic.
    
    Args:
        filename: Path to .par file
        
    Returns:
        Dictionary with parameter names as keys
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If file format is invalid
    """
    params: Dict[str, Union[float, List[float], List[int]]] = {}
    ltep_values: List[float] = []
    harm_values: List[int] = []
    lpap_triplets: List[List[float]] = []  # For tidal constituent triplets
    unmatched_periods: List[float] = []
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Skip empty lines
            line = line.strip()
            if not line:
                continue
                
            # Split into tag and values
            parts = line.split()
            if len(parts) < 2:
                continue
                
            tag = parts[0]
            values = parts[1:]
            
            try:
                if tag == "ltep":
                    # Long-term evolution parameters (LTE phases)
                    # Accumulate all ltep values into array
                    for val in values:
                        ltep_values.append(float(val))
                        
                elif tag == "harm":
                    # Harmonic indices (integers)
                    # Accumulate all harm values into array
                    for val in values:
                        harm_values.append(int(val))
                        
                else:
                    # Check if this is a numeric tag (period value for LPAP)
                    # Lines with numeric tags are tidal constituent triplets
                    try:
                        # Try parsing first value as float - if it's a number, 
                        # this is likely a period/amplitude/phase triplet
                        file_period = float(tag)
                        if len(values) == 2:
                            # This is a tidal constituent: period, amplitude, phase
                            # Tag is the period, values are amplitude and phase
                            amplitude = float(values[0])
                            phase = float(values[1])
                            
                            # Match to predefined LP period (within 1% tolerance)
                            matched_period = match_period(file_period)
                            if matched_period is not None:
                                # Use the predefined period, not the file period
                                lpap_triplets.append([matched_period, amplitude, phase])
                            else:
                                # No match found - log warning but include with file period
                                unmatched_periods.append(file_period)
                                lpap_triplets.append([file_period, amplitude, phase])
                                print(f"Warning: Line {line_num}: Period {file_period} "
                                      f"does not match any predefined LP period within 1%",
                                      file=sys.stderr)
                        elif len(values) == 1:
                            # Single value with numeric tag - treat as simple parameter
                            params[tag] = float(values[0])
                        else:
                            # Multiple values - store as array
                            params[tag] = [float(v) for v in values]
                    except ValueError:
                        # Not a numeric tag - standard parameter
                        if len(values) == 1:
                            params[tag] = float(values[0])
                        else:
                            # Multiple values - store as array
                            params[tag] = [float(v) for v in values]
                        
            except ValueError as e:
                print(f"Warning: Line {line_num}: Could not parse '{line}': {e}", 
                      file=sys.stderr)
                continue
    
    # Add accumulated arrays
    if ltep_values:
        params["ltep"] = ltep_values
    if harm_values:
        params["harm"] = harm_values
    if lpap_triplets:
        # Store LPAP as array of triplets (compatible with Ada JSON reader)
        params["lpap"] = lpap_triplets
    
    # Report unmatched periods as a warning
    if unmatched_periods:
        print(f"\nWARNING: {len(unmatched_periods)} periods did not match "
              f"predefined LP values:", file=sys.stderr)
        for p in unmatched_periods:
            print(f"  {p}", file=sys.stderr)
        print("These may not load correctly from JSON format.", file=sys.stderr)
        
    return params


def write_json(params: Dict, filename: str, indent: int = 2) -> None:
    """
    Write parameters to JSON file with maximum floating-point precision.
    
    Uses Python's repr() to preserve full IEEE 754 double precision.
    
    Args:
        params: Parameter dictionary
        filename: Output JSON file path
        indent: JSON indentation level (default: 2)
    """
    with open(filename, 'w') as f:
        # Encode with maximum float precision (all 17 significant digits)
        # This ensures exact reproduction of float values
        content = json.dumps(params, indent=indent, ensure_ascii=False)
        f.write(content)
        f.write('\n')  # Trailing newline


def convert_par_to_json(input_file: str, output_file: str = None) -> str:
    """
    Convert .par file to JSON format.
    
    Args:
        input_file: Input .par file path
        output_file: Output JSON file path (optional)
        
    Returns:
        Path to output file
        
    Raises:
        FileNotFoundError: If input file doesn't exist
    """
    input_path = Path(input_file)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Determine output filename
    if output_file is None:
        # Use .p extension (e.g., lt.exe.par -> lt.exe.p)
        output_path = input_path.with_suffix('.p')
    else:
        output_path = Path(output_file)
    
    # Parse .par file
    print(f"Reading {input_path}...", file=sys.stderr)
    params = parse_par_file(str(input_path))
    
    # Write JSON
    print(f"Writing {output_path}...", file=sys.stderr)
    write_json(params, str(output_path))
    
    # Report statistics
    scalar_count = sum(1 for v in params.values() if isinstance(v, (int, float)))
    array_count = sum(1 for v in params.values() if isinstance(v, list))
    
    print(f"Converted {scalar_count} scalar parameters, {array_count} arrays", 
          file=sys.stderr)
    
    if "ltep" in params:
        print(f"  - ltep: {len(params['ltep'])} values", file=sys.stderr)
    if "harm" in params:
        print(f"  - harm: {len(params['harm'])} values", file=sys.stderr)
    if "lpap" in params:
        print(f"  - lpap: {len(params['lpap'])} tidal constituents", file=sys.stderr)
    
    return str(output_path)


def main():
    """Main entry point for command-line usage."""
    if len(sys.argv) < 2:
        print("Usage: par_to_json.py <input.par> [output.json]", file=sys.stderr)
        print("", file=sys.stderr)
        print("Converts GEM-LTE .par parameter files to JSON format", file=sys.stderr)
        print("If output is not specified, uses input basename with .p extension", 
              file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        output_path = convert_par_to_json(input_file, output_file)
        print(f"Success: {output_path}", file=sys.stderr)
        sys.exit(0)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
