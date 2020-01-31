#!/usr/bin/python
import sys
import pandas as pd

if len(sys.argv) != 4:
    sys.stderr.write("Usage: %s <transmission_trees> <summary_stats> <alpha>\n" % sys.argv[0])
    sys.exit(1)

ttrees_filename = sys.argv[1]
summary_filename = sys.argv[2]
alpha = float(sys.argv[3])

#print(sys.argv)

df_summary = pd.read_csv(summary_filename, sep="\t")
df_summary = df_summary.sort_values(["unsampledLineages", "transmissions"])

cut_off = int(alpha * len(df_summary))
solIndices = set(df_summary["solIdx"][:cut_off+1])

with open(ttrees_filename) as f:
    f.readline()
    printedTrees = 0
    print(len(solIndices), "# trans trees")
    while printedTrees < len(solIndices):
        s = f.readline().rstrip("\n").split()
        numEdges = int(s[0])
        solIdx = int(s[-1])
        if solIdx in solIndices:
            printedTrees += 1
            print(" ".join(s))
            for i in range(numEdges):
                print(f.readline().rstrip("\n"))
        else:
            for i in range(numEdges):
                f.readline()
