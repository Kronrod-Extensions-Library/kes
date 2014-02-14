import sys
from numpy import array
from matplotlib.pyplot import *


if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: "+f)

F = open(f, "r")

found = False
matrix = []

for line in F.readlines():
    if not found and line.startswith(5*"="):
        found = True
        continue

    if found and line.startswith(5*"="):
        found = False
        continue

    if found:
        matrix.append(line.translate(None, "[]\n"))

F.close()

if not matrix:
    raise ValueError("No matrix data found!")

M = array([ map(int, line.split(" ")) for line in matrix if len(line) > 0 ])
maxn, maxnu = M.shape

print("Maximal values are: n="+str(maxn)+" and nu="+str(maxnu))

fig = figure()
imshow(M, interpolation="none", origin="upper", extent=[1,maxnu+1,maxn+1,1])
xlabel(r"$\nu$")
ylabel(r"$n$")
savefig("kronrod_extensions_flint_"+str(maxn)+"_"+str(maxnu)+".png")
