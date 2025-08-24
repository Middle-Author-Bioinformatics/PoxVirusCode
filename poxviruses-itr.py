import os
import re
from collections import defaultdict
from ArkTools import *

from Bio.Seq import Seq

def find_itr_length(genome_seq, min_len=10, max_len=100000, allowed_mismatches=10):
    """
    Identifies the length of inverted terminal repeats (ITRs) in a genome sequence.

    Parameters:
        genome_seq (str): The input DNA sequence.
        min_len (int): Minimum ITR length to consider.
        max_len (int): Maximum ITR length to search for.
        allowed_mismatches (int): Maximum number of mismatches allowed between ends.

    Returns:
        int: Length of the ITRs (0 if none found).
        str: ITR sequence from the left end.
    """
    genome = Seq(genome_seq.upper())
    max_check_len = min(max_len, len(genome) // 2)

    best_match_len = 0

    for itr_len in range(min_len, max_check_len + 1):
        left = genome[0:itr_len]
        right = genome[-itr_len:].reverse_complement()

        # Count mismatches
        mismatches = sum(1 for a, b in zip(left, right) if a != b)

        if mismatches <= allowed_mismatches:
            best_match_len = itr_len
        else:
            break

    return best_match_len

'''
signalDict = defaultdict(lambda: defaultdict(lambda: '-'))
signalp = open("/Users/agarber4/Downloads/signalP-pox.summary.csv")
for i in signalp:
    ls = i.rstrip().split(",")
    if ls[0] != "genome_assembly":
        assembly = ls[0].split(".ts")[0]
        signalDict[assembly]["secreted"] = ls[1]
        signalDict[assembly]["secreted with tm domain"] = ls[2]
        signalDict[assembly]["cytoplasmic"] = ls[3]
        signalDict[assembly]["cytoplasmic with tm domain"] = ls[4]

out = open("/Users/agarber4/Downloads/combined_stats.names.signalP.csv", "w")
summary = open("/Users/agarber4/Downloads/combined_stats.names.csv")
for i in summary:
    ls = i.rstrip().split(",")
    if ls[1] == "coding_density":
        out.write(i.rstrip() + ",secreted,secreted with tm domain,cytoplasmic,cytoplasmic with tm domain\n")
    else:
        sec = signalDict[ls[0]]["secreted"]
        secTM = signalDict[ls[0]]["secreted with tm domain"]
        cyto = signalDict[ls[0]]["cytoplasmic"]
        cytoTM = signalDict[ls[0]]["cytoplasmic with tm domain"]
        out.write(i.rstrip() + "," + str(sec) + "," + str(secTM) + "," + str(cyto) + "," + str(cytoTM) + "\n")
out.close()
'''

contigDict = defaultdict(lambda: '-')
contigs = os.listdir("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/contigs")
for i in contigs:
    if re.search(r'fa', i):
        file = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/contigs/%s" % i)
        file = fasta(file)
        assembly = i.split(".f")[0]
        if len(file.keys()) == 1 and re.findall(r'complete', list(file.keys())[0]):
            for j in file.keys():
                contigDict[assembly] = file[j]


totalDict = defaultdict(list)
parDict = defaultdict(lambda: defaultdict(list))
dupDict = defaultdict(lambda: defaultdict(list))
blastpDir = os.listdir("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/blastp")
for i in blastpDir:
    if re.search(r'tsv', i):
        file = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/blastp/%s" % i)
        for j in file:
            ls = j.rstrip().split("\t")
            if ls[0] != ls[1]:
                pair = sorted([ls[0], ls[1]])
                if pair[0] not in dupDict[i]:
                    if float(ls[2]) == 100:
                        dupDict[i][pair[0]].append(pair[1])
                    else:
                        parDict[i][pair[0]].append(pair[1])
            else:
                totalDict[i].append(ls[0])

out = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/ITR.duplicates.paralogs.singletons.csv", "w")
out.write("genome,ITR_length,duplicate_clusters,total_duplicates,paralog_clusters,total_paralogs,total_genes\n")
for i in dupDict.keys():
    assembly = i.split(".blas")[0]
    if assembly in contigDict.keys():

        seq = (contigDict[assembly])
        print(i)
        itr = (find_itr_length(seq))
        pars = (len(parDict[i]))
        dups = (len(dupDict[i]))
        tots = (len(totalDict[i]))
        totalDups = 0
        totalPars = 0
        for j in dupDict[i]:
            totalDups += len(dupDict[i][j])

        for j in parDict[i]:
            totalPars += len(parDict[i][j])
        out.write(assembly + "," + str(itr) + "," + str(dups) + "," + str(totalDups) + "," + str(pars) + "," + str(totalPars) + "," + str(tots) + "\n")

out.close()
















