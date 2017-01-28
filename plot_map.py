import sys
import re
from numpy import array
from matplotlib.pyplot import figure, imshow, grid, xlabel, ylabel, tight_layout, savefig, cm


if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: {}".format(f))
else:
    raise ValueError("No file given!")

found = False
matrix = []

with open(f, "r") as F:
    for line in F.readlines():
        if not found and line.startswith(5 * "="):
            found = True
            continue

        if found and line.startswith(5 * "="):
            found = False
            continue

        if found:
            matrix.append(re.sub('[[\]\n]', '', line))

if not matrix:
    raise ValueError("No matrix data found!")

M = array([[int(c) for c in line.split(" ")] for line in matrix if len(line) > 0])
maxn, maxp = M.shape

print("Maximal values are: n={} and p={}".format(maxn, maxp))

figure()
imshow(M, interpolation="none", origin="upper", extent=[1, maxp + 1, maxn + 1, 1], cmap=cm.binary)
grid(True)
xlabel(r"$p$")
ylabel(r"$n$")
tight_layout()
savefig(f[:-4] + ".png")
