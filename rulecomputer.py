import sys
import concurrent.futures
import subprocess

from rulelistfile import get_rulelist


def compute_rules(rulelist, call, max_workers=8):

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:

        F = []

        for rule in rulelist:

            def compute_rule(rule):

                def doit():
                    ruledatafile = "rule_{}.txt".format("_".join(map(str, rule)))
                    with open(ruledatafile, 'a') as stdout:
                        subprocess.call(call(rule), shell=False, stdout=stdout)

                return doit

            F.append(executor.submit(compute_rule(rule)))

        for future in concurrent.futures.as_completed(F):
            future.result()


if __name__ == '__main__':

    call = lambda rule: ['./kes', '-cn', '-cw', *map(str, rule)]

    if not len(sys.argv[1:]) >= 1:
        raise IOError('No input files given to read rules from.')

    # Merge all rules before calling computep
    rules = [rule for listfile in sys.argv[1:] for rule in get_rulelist(listfile)]

    # Compute all rules in parallel
    compute_rules(rules, call)
