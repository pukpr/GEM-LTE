#!/usr/bin/env python3
"""Verify JSON file structure"""
import json
import sys

with open('lt.exe.p', 'r') as f:
    data = json.load(f)

print("JSON structure:")
print("=" * 60)
print(f"Number of keys: {len(data)}")
print(f"Keys: {list(data.keys())}")
print()

# Check all expected fields
expected_scalars = ['offs', 'bg', 'impA', 'impB', 'impC', 'delA', 'delB', 
                    'asym', 'ann1', 'ann2', 'sem1', 'sem2', 'year', 'IR', 
                    'ma', 'mp', 'shfT', 'init']
print(f"Scalar fields ({len(expected_scalars)} expected):")
for field in expected_scalars:
    if field in data:
        print(f"  ✓ {field}: {data[field]}")
    else:
        print(f"  ✗ {field}: MISSING!")

print()
print(f"Array fields:")
print(f"  ltep: {len(data.get('ltep', []))} values")
print(f"  harm: {len(data.get('harm', []))} values") 
print(f"  lpap: {len(data.get('lpap', []))} triplets")

print()
print("First 3 LPAP triplets:")
for i, trip in enumerate(data.get('lpap', [])[:3]):
    print(f"  {i+1}. period={trip[0]:.11f}, amp={trip[1]:.11f}, phase={trip[2]:.11f}")

print()
print("LPAP array structure check:")
if 'lpap' in data:
    lpap = data['lpap']
    print(f"  Type: {type(lpap)}")
    print(f"  Length: {len(lpap)}")
    if len(lpap) > 0:
        print(f"  First element type: {type(lpap[0])}")
        if isinstance(lpap[0], list):
            print(f"  First element length: {len(lpap[0])}")
