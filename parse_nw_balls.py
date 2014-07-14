import sys
import os
import re
from numpy import array, savez

if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: "+f)
else:
    raise ValueError("No file given!")

rblock = False
wblock = False
allroots = []
allweights = []
roots = []
weights = []

with open(f, "r") as F:
    for line in F.readlines():
        if line.startswith("The nodes are"):
            rblock = True
            wblock = False
            roots = []
            allroots.append(roots)
            continue
        elif line.startswith("The weights are"):
            wblock = True
            rblock = False
            weights = []
            allweights.append(weights)
            continue
        elif (rblock or wblock) and not line.startswith("| "):
            rblock = False
            wblock = False
            continue
        elif rblock:
            N = map(float, re.findall("[-+]?\d+\.?\d*[Ee]?[+-]?\d*", line))
            roots.append(N[0] + 1.0j*N[1])
        elif wblock:
            N = map(float, re.findall("[-+]?\d+\.?\d*[Ee]?[+-]?\d*", line))
            weights.append(N[0] + 1.0j*N[1])

if not allroots or not allweights:
    raise ValueError("No suitable data found!")

allroots = map(array, allroots)
allweights = map(array, allweights)

b = os.path.basename(f).rsplit(".",1)[0]

for i, r in enumerate(allroots):
    with open(b+"_nodes_"+str(i)+".dat", "w") as f:
        savez(f, r)

for i, w in enumerate(allweights):
    with open(b+"_weights_"+str(i)+".dat", "w") as f:
        savez(f, w)
