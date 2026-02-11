import sympy as sp
import csv

def read_cross_harmonics(filename):
    harmonics = []
    with open(filename, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            i, j, k, amp, phase = row
            harmonics.append((int(i), int(j), int(k), float(amp), float(phase)))
    return harmonics

def only_A_C_crossharmonic(expr, A, B, C):
    """
    Rewrites expr so that only cross-harmonics of form cos(iA + kC) or sin(iA + kC) remain.
    Any other cross-harmonic involving B is expanded into products/sums of pure harmonics.
    """
    # Recursively process all trig functions in expr
    expr = sp.expand_trig(expr)
    expr = expr.replace(
        lambda e: (
            (isinstance(e, sp.sin) or isinstance(e, sp.cos))
            and e.args
            and e.args[0].has(B)
        ),
        lambda e: _project_out_B(e, A, B, C)
    )
    return sp.simplify(expr)

def _project_out_B(trig_expr, A, B, C):
    """
    For a trig_expr of form cos(iA + jB + kC), expand so that B is not in a cross-harmonic argument.
    """
    x = trig_expr.args[0]
    # Try to decompose x into iA + jB + kC + phase
    x = sp.expand(x)
    coeffs = {A: 0, B: 0, C: 0}
    phase = 0
    if isinstance(x, sp.Add):
        for term in x.args:
            if term.has(A):
                coeffs[A] += int(term.as_coeff_exponent(A)[0])
            elif term.has(B):
                coeffs[B] += int(term.as_coeff_exponent(B)[0])
            elif term.has(C):
                coeffs[C] += int(term.as_coeff_exponent(C)[0])
            else:
                phase += term
    else:
        # Single term
        if x.has(A):
            coeffs[A] += int(x.as_coeff_exponent(A)[0])
        elif x.has(B):
            coeffs[B] += int(x.as_coeff_exponent(B)[0])
        elif x.has(C):
            coeffs[C] += int(x.as_coeff_exponent(C)[0])
        else:
            phase += x

    i, j, k = coeffs[A], coeffs[B], coeffs[C]
    # If j == 0, leave unchanged (no B in argument)
    if j == 0:
        return trig_expr.func(i*A + k*C + phase)
    # Otherwise, expand using sum-of-angles with B isolated
    arg1 = i*A + k*C + phase
    arg2 = j*B
    if isinstance(trig_expr, sp.cos):
        # cos(arg1 + arg2) = cos(arg1)*cos(arg2) - sin(arg1)*sin(arg2)
        return sp.expand_trig(sp.cos(arg1)*sp.cos(arg2) - sp.sin(arg1)*sp.sin(arg2))
    elif isinstance(trig_expr, sp.sin):
        # sin(arg1 + arg2) = sin(arg1)*cos(arg2) + cos(arg1)*sin(arg2)
        return sp.expand_trig(sp.sin(arg1)*sp.cos(arg2) + sp.cos(arg1)*sp.sin(arg2))
    else:
        return trig_expr  # If not sin/cos

def canonicalize_expr(expr, A, B, C):
    """
    Canonicalize expr to allow only pure harmonics and A/C cross-harmonics.
    """
    # Expand and simplify
    expr = sp.expand_trig(expr)
    expr = only_A_C_crossharmonic(expr, A, B, C)
    expr = sp.expand_trig(expr)
    expr = sp.simplify(expr)
    expr = sp.collect(expr, [sp.sin(A), sp.cos(A), sp.sin(B), sp.cos(B), sp.sin(C), sp.cos(C)])
    return expr

def main(input_file):
    A, B, C = sp.symbols('A B C', real=True)
    harmonics = read_cross_harmonics(input_file)
    total_expr = 0
    for i, j, k, amplitude, phase in harmonics:
        term = amplitude * sp.cos(i*A + j*B + k*C + phase)
        total_expr += term
        print("term")

    canonical_expr = canonicalize_expr(total_expr, A, B, C)
    print("Canonical symbolic expression (only fundamentals/harmonics, A/C cross-harmonics allowed):")
    print(canonical_expr)
    print("\nPretty print:")
    # sp.pprint(canonical_expr, use_unicode=True)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python canonicalize_fundamental_harmonics.py input.csv")
        exit(1)
    main(sys.argv[1])