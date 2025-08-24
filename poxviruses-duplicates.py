import os
import re
import statistics
from collections import defaultdict

'''
lengthDict = defaultdict(lambda: '-')
lengths = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/lengths.txt")
for i in lengths:
    ls = i.rstrip().split(" ")
    lengthDict[ls[0].split(".")[0]] = float(ls[1])

poxDict = defaultdict(lambda: defaultdict(list))
poxDir = os.listdir("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/blast")
for i in poxDir:
    if re.findall(r'blast', i):
        accession = (i.split(".")[0])
        if accession in lengthDict.keys():
            length = (lengthDict[accession])
            file = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/blast/" + i)
            for j in file:
                ls = j.rstrip().split("\t")
                qmean = statistics.mean([int(ls[6]), int(ls[7])])
                smean = statistics.mean([int(ls[8]), int(ls[9])])
                diff = qmean - smean
                if diff != 0 and int(ls[3]) > 1000:
                    prop = diff / length
                    if prop < -0.5 or prop > 0.5:
                        poxDict[accession][ls[3]].append(float(ls[2]))

out = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/ITR.csv", "w")
out.write("accession,genome_length,ITR_size,mean_ANI\n")
for i in poxDict.keys():
    genomeLength = (lengthDict[i])
    for j in poxDict[i]:
        totalLength = 0
        IDs = []
        if len(poxDict[i][j]) == 2:
            length = int(j)
            totalLength += length
            ID = float(poxDict[i][j][0])
            IDs.append(str(ID))
            IDs = IDs*length
        IDs = [float(x) for x in IDs]
        try:
            out.write(i + "," + str(genomeLength) + "," + str(totalLength) + "," + str(statistics.mean(IDs)) + "\n")
        except statistics.StatisticsError:
            pass

nameDict = defaultdict(lambda: defaultdict(lambda: 'na'))
ncbi = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/ARB_review/ncbi_assembly_info.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        nameDict[ls[0]]["name"] = ls[7]
        nameDict[ls[0]]["strain"] = ls[8]


summary = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combined_stats.tsv")
out = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combined_stats.names.tsv", "w")
# out.write("assembly\tcoding_density\tcoverage\tcontigs\tcontig_length\tITR_size\tmean_ANI\tgenome_length\tITR_size\tmean_ANI\tname\tstrain\n")
for i in summary:
    ls = i.rstrip().split("\t")
    out.write(i.rstrip() + "\t" + str(nameDict[ls[0]]["name"]) + "\t" + str(nameDict[ls[0]]["strain"]) + "\n")
out.close()
'''

itrDict = defaultdict(lambda: defaultdict(lambda: '-'))
itrs = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/ITR.duplicates.paralogs.singletons.csv")
for i in itrs:
    ls = i.rstrip().split(",")
    itrDict[ls[0]]["itr_length"] = ls[1]
    itrDict[ls[0]]["duplicate_clusters"] = ls[2]
    itrDict[ls[0]]["total_duplicates"] = ls[3]
    itrDict[ls[0]]["paralog_clusters"] = ls[4]
    itrDict[ls[0]]["total_paralogs"] = ls[5]
    itrDict[ls[0]]["total_genes"] = ls[6]

dateDict = defaultdict(lambda: '-')
dates = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/ncbi_assembly_info.pox.tsv")
for i in dates:
    ls = i.rstrip().split("\t")
    dateDict[ls[0]] = ls[14]

out = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/combined_stats.names.signalP.dates.itrs.duplicates.paralogs.singletons.csv", "w")
combined = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/combined_stats.names.signalP.csv")
for i in combined:
    ls = i.rstrip().split(",")
    if ls[1] == "coding_density":
        out.write(i.rstrip() + ",itr_length,duplicate_clusters,total_duplicates,paralog_clusters,total_paralogs,total_genes,date\n")
    else:
        if ls[0] in itrDict.keys():
            out.write(i.rstrip() + "," + str(itrDict[ls[0]]["itr_length"]) + "," + str(itrDict[ls[0]]["duplicate_clusters"]) + "," + str(itrDict[ls[0]]["total_duplicates"]) + "," + str(itrDict[ls[0]]["paralog_clusters"]) + "," + str(itrDict[ls[0]]["total_paralogs"]) + "," + str(itrDict[ls[0]]["total_genes"]) + "," + str(dateDict[ls[0]]) + "\n")
        else:
            out.write(i.rstrip() + ",-,-,-,-,-,-,-\n")
out.close()


curated = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/File_1.csv")
out = open("/Users/agarber4/Desktop/Ongoing_Research_Projects/poxviruses/combine/File_1.allInfo.csv", "w")
for i in curated:
    ls = i.rstrip().split(",")
    if ls[1] == "coding_density":
        out.write(i.rstrip() + ",itr_length,duplicate_clusters,total_duplicates,paralog_clusters,total_paralogs,total_genes,date\n")
    else:
        if ls[0] in itrDict.keys():
            out.write(i.rstrip() + "," + str(itrDict[ls[0]]["itr_length"]) + "," + str(itrDict[ls[0]]["duplicate_clusters"]) + "," + str(itrDict[ls[0]]["total_duplicates"]) + "," + str(itrDict[ls[0]]["paralog_clusters"]) + "," + str(itrDict[ls[0]]["total_paralogs"]) + "," + str(itrDict[ls[0]]["total_genes"]) + "," + str(dateDict[ls[0]]) + "\n")
        else:
            out.write(i.rstrip() + ",-,-,-,-,-,-,-\n")
out.close()




