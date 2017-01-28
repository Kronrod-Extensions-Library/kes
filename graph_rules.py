import sys
from functools import reduce


if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: {}".format(f))
else:
    raise ValueError("No file given!")

ruletree = {}

with open(f, "r") as F:
    for line in F.readlines():
        if line.startswith("Search for recursive extensions of:"):
            n = line[37:-1]
            print("n: {}".format(n))
        if line.startswith("Maximal allowed extension order p: "):
            p = line[35:-1]
            print("Max. p: {}".format(p))
        if line.startswith("Maximal allowed recursion depth: "):
            rec = line[33:-1]
            print("Max. depth: {}".format(rec))

        if line.startswith("RULE:"):
            print(line[:-1])
            L = line.split()
            data = L[2:(4+int(L[1]))]
            t = ruletree
            for digit in data:
                if digit in t:
                    t = t[digit]
                else:
                    t[digit] = {}

def print_rule(rule, rd):
    for k, v in rd.items():
        if rule:
            ruleoldstr = "\"Q[" + reduce(lambda x,y: str(x)+","+str(y), rule) + "]\""
        else:
            ruleoldstr = "\"\""
        rule.append(k)
        rulenewstr = "\"Q[" + reduce(lambda x,y: str(x)+","+str(y), rule) + "]\""
        graph.write(rulenewstr + " [shape=box];\n")
        graph.write(ruleoldstr + " -> " + rulenewstr + ";\n")
        if v:
            print_rule(rule, v)
        rule.pop()

gf = "graph_{}_{}_{}.dot".format(n, p, rec)
print("Writing graph to file: {}".format(gf))

with open(gf, "w") as graph:
    graph.write("digraph G {\n")
    graph.write("  rankdir=LR;")

    rule = []
    t = ruletree
    print_rule(rule, ruletree)

    graph.write("}\n")
