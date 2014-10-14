import sys
import re
from numpy import array, savez

if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: "+f)
else:
    raise ValueError("No file given!")

rblock = False
wblock = False
roots = []
weights = []

with open(f, "r") as F:
    for line in F.readlines():
        if line.startswith("DIMENSION"):
            dimension = int(re.findall("\d+", line)[0])
        elif line.startswith("LEVEL"):
            level = int(re.findall("\d+", line)[0])
        elif line.startswith("NODES"):
            rblock = True
            continue
        elif line.startswith("WEIGHTS"):
            wblock = True
            continue
        elif (rblock or wblock) and not line.startswith("| "):
            rblock = False
            wblock = False
            continue
        elif rblock:
            parts = line.split(",")
            parts = map(lambda t: t.strip(), parts)
            parts = [p for p in parts if len(p) > 0]
            matches = map(lambda t: re.findall("[-+]?\d+\.?\d*[Ee]?[+-]?\d*", t), parts)
            numbers = [map(float, u) for u in matches]
            roots.append([n[0]+1j*n[1] for n in numbers])
        elif wblock:
            N = map(float, re.findall("[-+]?\d+\.?\d*[Ee]?[+-]?\d*", line))
            weights.append(N[0] + 1.0j*N[1])

if not roots or not weights:
    raise ValueError("No suitable data found!")

print("Dimension: %d" % dimension)
print("Level: %d" % level)

nr = len(roots)
nw = len(weights)
roots = array(roots)
weights = array(weights)
assert nr == nw

with open("gk_nodes_dimension_%d_level_%d.dat" % (dimension, level), "w") as f:
    savez(f, roots)

with open("gk_weights_dimension_%d_level_%d.dat" % (dimension, level), "w") as f:
    savez(f, weights)
