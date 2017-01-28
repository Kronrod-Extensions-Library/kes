import sys
import re
from numpy import array, savez

if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: {}".format(f))
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
            parts = [t.strip() for t in parts]
            parts = [p for p in parts if len(p) > 0]
            matches = [re.findall("[-+]?\d+\.?\d*[Ee]?[+-]?\d*", t) for t in parts]
            numbers = [[float(ui) for ui in u] for u in matches]
            roots.append([n[0]+1j*n[1] for n in numbers])
        elif wblock:
            N = [float(c) for c in re.findall("[-+]?\d+\.?\d*[Ee]?[+-]?\d*", line)]
            weights.append(N[0] + 1.0j*N[1])

if not roots or not weights:
    raise ValueError("No suitable data found!")

print("Dimension: {:d}".format(dimension))
print("Level: {:d}".format(level))

nr = len(roots)
nw = len(weights)
roots = array(roots)
weights = array(weights)
assert nr == nw

filename = "gk_nodes_dimension_{:d}_level_{:d}".format(dimension, level)
savez(filename, roots)

filename = "gk_weights_dimension_{:d}_level_{:d}".format(dimension, level)
savez(filename, weights)
