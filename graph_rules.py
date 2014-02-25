import sys


if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: "+f)

F = open(f, "r")

ruletree = {}

for line in F.readlines():
    if line.startswith("Search for recursive extensions of: P"):
        n = line[37:-1]
        print("n: "+str(n))
    if line.startswith("Maximal allowed extension order p: "):
        p = line[35:-1]
        print("Max. p: "+str(p))
    if line.startswith("Maximal allowed recursion depth: "):
        rec = line[33:-1]
        print("Max. depth: "+str(rec))

    if line.startswith("RULE:"):
        print(line[:-1])
        L = line.split()
        data = L[4:(4+int(L[1]))]
        t = ruletree
        for digit in data:
            if t.has_key(digit):
                t = t[digit]
            else:
                t[digit] = {}

F.close()

# Color rules not found with old code
even = lambda x: int(x) % 2 == 0
fc = lambda l: {True:"black", False:"red"}[all(map(even, l))]

def print_rule(rule, rd):
    for k, v in rd.iteritems():
        if rule:
            ruleoldstr = "\"Q[" + reduce(lambda x,y: str(x)+","+str(y), rule) + "]\""
        else:
            ruleoldstr = "\"\""
        rule.append(k)
        rulenewstr = "\"Q[" + reduce(lambda x,y: str(x)+","+str(y), rule) + "]\""
        graph.write(rulenewstr + " [shape=box,color="+fc(rule[1:])+"];\n")
        graph.write(ruleoldstr + " -> " + rulenewstr + ";\n")
        if v:
            print_rule(rule, v)
        rule.pop()

gf = "graph_"+str(n)+"_"+str(p)+"_"+str(rec)+".dot"
print("Writing graph to file: "+gf)

graph = file(gf, "w")
graph.write("digraph G {\n")
graph.write("  rankdir=LR;")

rule = []
t = ruletree
print_rule(rule, ruletree)

graph.write("}\n")
graph.close()
