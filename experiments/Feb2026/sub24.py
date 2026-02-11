import sys

delta = 1.0/24.0

fname = sys.argv[1]

lines = open(fname).readlines()
out = open(fname, "w")

for line in lines:
    f = line.split()
    f[0] = "%.7F" % (float(f[0]) - delta)
    out.write("    ".join(f) + "\n")

out.close()
