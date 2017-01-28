import sys
import os
import re
from numpy import array, savez, savetxt, real, imag

if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: {}".format(f))
else:
    raise IOError("No file given!")

rblock = False
wblock = False
allroots = []
allweights = []
allrootserr = []
allweightserr = []

with open(f, "r") as F:
    for line in F.readlines():
        if line.startswith("The nodes are"):
            rblock = True
            wblock = False
            roots = []
            rootserr = []
            allroots.append(roots)
            allrootserr.append(rootserr)
            continue
        elif line.startswith("The weights are"):
            wblock = True
            rblock = False
            weights = []
            weightserr = []
            allweights.append(weights)
            allweightserr.append(weightserr)
            continue
        elif (rblock or wblock) and not line.startswith("| "):
            rblock = False
            wblock = False
            continue
        elif rblock:
            N = [float(s.replace(' ', '')) for s in re.findall("[-+]?[ ]?\d+\.?\d*[Ee]?[+-]?\d*", line)]
            roots.append(N[0] + 1.0j * N[1])
            rootserr.append(N[2] + 1.0j * N[3])
        elif wblock:
            W = [float(s.replace(' ', '')) for s in re.findall("[-+]?[ ]?\d+\.?\d*[Ee]?[+-]?\d*", line)]
            weights.append(W[0] + 1.0j * W[1])
            weightserr.append(W[2] + 1.0j * W[3])

if not allroots or not allweights:
    raise ValueError("No suitable data found!")

allroots = map(array, allroots)
allweights = map(array, allweights)
allrootserr = map(array, allrootserr)
allweightserr = map(array, allweightserr)

b = os.path.basename(f).rsplit(".", 1)[0]

for i, r in enumerate(allroots):
    filename = b + "_nodes_" + str(i)
    savez(filename, r)

for i, w in enumerate(allweights):
    filename = b + "_weights_" + str(i)
    savez(filename, w)

for i, (r, w, rerr, werr) in enumerate(zip(allroots, allweights, allrootserr, allweightserr)):
    with open(b + "_qr_" + str(i) + ".csv", "w") as f:
        savetxt(f, array([real(r), imag(r),
                          real(w), imag(w),
                          real(rerr), imag(rerr),
                          real(werr), imag(werr)]), delimiter=", ")
