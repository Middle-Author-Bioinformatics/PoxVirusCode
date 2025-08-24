from collections import defaultdict

nameDict = defaultdict(lambda: '-')

names = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combined_stats.names.tsv")
for i in names:
    ls = i.rstrip().split("\t")
    nameDict[ls[0].split(".")[0]] = ls[11].split(" ")[0]

out = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/codon_usages/codon_usages.names.tsv", "w")
pox = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/codon_usages/codon_usages.tsv")
for i in pox:
    ls = i.rstrip().split("\t")
    if ls[1] != "GC":
        if ls[0] in nameDict.keys():
            name = nameDict[ls[0]]
            if name == "Yaba-like":
                name = "Yaba"
            out.write(i.rstrip() + "\t" + str(name) + "\n")
    else:
        out.write(i.rstrip() + "\tviral_taxa\n")
out.close()