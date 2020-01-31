#!/usr/bin/python
import sys
import copy

def parse_dimacs_vars(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            s = line.rstrip("\n").split()
            if len(s) > 3 and s[3] == "(v,":
                d[int(s[0])] = map(int, s[1:3])
    return d

def parse_unigen(filename, max_var):
    d = {}
    with open(filename) as f:
        for line in f:
            s = line.lstrip(" ").lstrip("v").rstrip("\n").split()
            if len(s) == 0: continue
            count = int(s[-1].split(":")[-1])
            assignment = frozenset(map(int, s[:-1])[:max_var])

            if assignment not in d:
                d[assignment] = count
            else:
                d[assignment] += count

        return d

def write_ptree(in_ptree, var_mapping, sol, prefix, idx, count):
    out_ptree = copy.deepcopy(in_ptree)

    for lit in sol:
        if lit > 0:
            v,s = var_mapping[lit]
            out_ptree[v-1][-1] = str(s+1)

    with open("%sidx%d_count%d.out" % (prefix, idx, count), "w") as f:
        for line in out_ptree:
            f.write("\t".join(line) + "\n")

def parse_ptree(filename):
    res = []
    with open(filename) as f:
        for line in f:
            res.append(line.rstrip("\n").split())
        return res

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: %s <ptree> <dimacs_vars> <unigen_sol> <prefix>\n" % sys.argv[0])
        sys.exit(1)

    var_mapping = parse_dimacs_vars(sys.argv[2])
    max_var = max(var_mapping.keys())

    #print var_mapping
    #print max_var

    in_ptree = parse_ptree(sys.argv[1])
    solutions = parse_unigen(sys.argv[3], max_var)
    prefix = sys.argv[4]
    
    for idx, sol in enumerate(solutions):
        write_ptree(in_ptree, var_mapping, sol, prefix, idx, solutions[sol])
