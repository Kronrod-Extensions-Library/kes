import os
import os.path as path
import re
import itertools


def get_rulelist(file):
    rules = []

    with open(file, 'r') as f:
        for line in f.readlines():
            if line.startswith('RULE:'):
                rule = re.split('[\s]+', line.strip())[2:]
                rules.append(tuple(map(int, rule)))

    return rules


def collect_all_rules(rulelistspath):

    allrules = {}

    for file in os.listdir(rulelistspath):
        print(file)
        if file.startswith('rules_'):
            m = re.match('rules_n(.*)_maxp(.*)_maxrec(.*).txt', file)
            datum = tuple(map(int, (m.group(1), m.group(2), m.group(3))))
            allrules[datum] = get_rulelist(path.join(rulelistspath, file))

    return allrules


def merge_rulelists(allrules):

    mergedrules = {}

    key = lambda x: x[0]
    for g, si in itertools.groupby(sorted(allrules.keys(), key=key), key):
        pmax = 0
        rmax = 0
        rules = []

        for i in si:
            n, p, r = i
            pmax = max(pmax, p)
            rmax = max(rmax, r)
            rules += allrules[i]

        mergedrules[(g, pmax, rmax)] = sorted(list(set(rules)))

    return mergedrules
