#!/usr/bin/env python3
"""
json_to_par.py - Convert GEM-LTE JSON files back to .par format

Reverses the par_to_json.py conversion to verify data completeness.
If this round-trip works correctly, it proves the JSON format contains
all necessary information and any loading issues are in the Ada JSON reader.

Usage:
    python json_to_par.py <input.json> [output.par]
    
If output filename is not specified, uses input basename with .par extension.
"""

import sys
import json
from pathlib import Path
from typing import Dict, List, Union


def write_par_file(params: Dict, filename: str) -> None:
    """
    Write parameters to .par file in tagged-value format.
    
    Args:
        params: Parameter dictionary from JSON
        filename: Output .par file path
    """
    with open(filename, 'w') as f:
        # Write scalar parameters
        scalars = ['offs', 'bg', 'impA', 'impB', 'impC', 'delA', 'delB',
                   'asym', 'ann1', 'ann2', 'sem1', 'sem2', 'year', 'IR',
                   'ma', 'mp', 'shfT', 'init']
        
        for key in scalars:
            if key in params:
                value = params[key]
                # Match original .par formatting: key with spaces, then value
                f.write(f"{key:4s} {value:17.11f}\n")
        
        # Write LPAP triplets (period, amplitude, phase)
        # Period is used as the tag
        if 'lpap' in params:
            for triplet in params['lpap']:
                period, amplitude, phase = triplet
                # Format: period as tag (14.11f), then amplitude and phase
                f.write(f"{period:17.11f} {amplitude:17.11f} {phase:17.11f}\n")
        
        # Write ltep array (long-term evolution parameters)
        if 'ltep' in params:
            for value in params['ltep']:
                f.write(f"ltep {value:17.11f}\n")
        
        # Write harm array (harmonic indices)
        if 'harm' in params:
            for value in params['harm']:
                f.write(f"harm {value:3d}\n")


def convert_json_to_par(input_file: str, output_file: str = None) -> str:
    """
    Convert JSON file to .par format.
    
    Args:
        input_file: Input JSON file path
        output_file: Output .par file path (optional)
        
    Returns:
        Output file path
    """
    input_path = Path(input_file)
    
    if output_file is None:
        # Default output: replace .p or .json with .par
        if input_path.suffix == '.p':
            output_file = str(input_path.with_suffix('.par'))
        elif input_path.suffix == '.json':
            # Remove .json suffix and add .par
            output_file = str(input_path.with_suffix('')) + '.par'
        else:
            output_file = str(input_path) + '.par'
    
    print(f"Reading {input_file}...", file=sys.stderr)
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    print(f"Writing {output_file}...", file=sys.stderr)
    write_par_file(params, output_file)
    
    # Report statistics
    scalar_count = sum(1 for k in params.keys() 
                      if k not in ['ltep', 'harm', 'lpap'])
    ltep_count = len(params.get('ltep', []))
    harm_count = len(params.get('harm', []))
    lpap_count = len(params.get('lpap', []))
    
    print(f"Converted {scalar_count} scalar parameters", file=sys.stderr)
    if ltep_count > 0:
        print(f"  - ltep: {ltep_count} values", file=sys.stderr)
    if harm_count > 0:
        print(f"  - harm: {harm_count} values", file=sys.stderr)
    if lpap_count > 0:
        print(f"  - lpap: {lpap_count} tidal constituents", file=sys.stderr)
    print(f"Success: {output_file}", file=sys.stderr)
    
    return output_file


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python json_to_par.py <input.json> [output.par]", 
              file=sys.stderr)
        print("", file=sys.stderr)
        print("Convert GEM-LTE JSON parameter file to .par format", 
              file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        convert_json_to_par(input_file, output_file)
    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {input_file}: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
