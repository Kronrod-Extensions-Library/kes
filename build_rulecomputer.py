import sys

ruletype = "hermite"

if len(sys.argv) == 2:
    f = sys.argv[1]
    print("Reading values from file: "+f)
else:
    raise ValueError("No file given!")

basecommand = "./kes -cn -cw "
commands = []

with open(f, "r") as F:
    for line in F.readlines():
        l = line[:-2]
        command = basecommand + l + " > rule_" + ruletype + "_" + l.replace(" ", "_") + ".dat\n"
        commands.append(command)

with open("compute_all_rules.sh", "w") as F:
    F.writelines(commands)
