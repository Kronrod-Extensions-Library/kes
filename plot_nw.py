import sys
import re
from numpy import array, real, imag
from matplotlib.pyplot import figure, plot, semilogy, grid, xlim, ylim, xlabel, ylabel, savefig


if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: {}".format(f))
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
        if line.startswith("Extension levels are"):
            L = line.split(": ")[1][:-1]
        elif line.startswith("The nodes are"):
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
            N = [float(s.replace(' ', '')) for s in re.findall("[-+]?[ ]?\d+\.?\d*[Ee]?[+-]?\d*", line)]
            roots.append(N[0] + 1.0j * N[1])
        elif wblock:
            W = [float(s.replace(' ', '')) for s in re.findall("[-+]?[ ]?\d+\.?\d*[Ee]?[+-]?\d*", line)]
            weights.append(W[0] + 1.0j * W[1])

if not allroots or not allweights:
    raise ValueError("No suitable data found!")

L = L.replace(" ", "_")
allroots = map(array, allroots)
allweights = map(array, allweights)


fig = figure()
for roots in allroots:
    plot(real(roots), imag(roots), "o")
    xlim(min(real(roots).min() - 1, -1), max(real(roots).max() + 1, 1))
    ylim(min(imag(roots).min(), -1), max(imag(roots).max(), 1))
grid(True)
xlabel(r"$\Re \gamma_i$")
ylabel(r"$\Im \gamma_i$")
savefig(f[:-4] + "_nodes_weights_cp_" + L + ".svg")


fig = figure()
for roots, weights in zip(allroots, allweights):
    plot(real(roots), real(weights), "o")
    xlim(min(real(roots).min() - 1, -1), max(real(roots).max() + 1, 1))
grid(True)
xlabel(r"$\gamma_i$")
ylabel(r"$\omega_i$")
savefig(f[:-4] + "_nodes_weights_" + L + ".svg")


fig = figure()
for roots, weights in zip(allroots, allweights):
    semilogy(real(roots), abs(real(weights)), "o")
    xlim(min(real(roots).min() - 1, -1), max(real(roots).max() + 1, 1))
grid(True)
xlabel(r"$\gamma_i$")
ylabel(r"$\omega_i$")
savefig(f[:-4] + "_nodes_weights_log_" + L + ".svg")
