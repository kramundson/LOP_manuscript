#!/usr/bin/env/ python3
# indel-allelic-ratio-3.py
# Kirk Amundson
# 16 April 2019

from optparse import OptionParser
from collections import defaultdict

"""
Program computes ratio of A vs B allele across user-provided intervals
Two files are needed:
    1. Alleles file, which reports allele-specific read depth for each sample at each
       polymorphic locus of interest
    2. Lesions file, which reports sample-specific intervals to compute allele ratio over.
    
Alleles file follows format of CallAllelesAB.py output. In tab-delimite text:
Chrom   Pos Ref allele_A    allele_B    Snptype-sample1 SNP1-sample1    SNP2-sample1 TotalCov-sample1   %A-sample1  CovA-sample1

Note, Snptype thru CovA then repeat for next individual in the dataset

Lesion file is also tab-delimited text:

chrom   lesion_start    lesion_end  sample_with_lesion  lesion_size

See included example files for testing
"""

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-a", "--allelecalls", dest="a", help="allele calls file")
parser.add_option("-l", "--lesions", dest="l", help="lesion file")
parser.add_option("-o", "--out", dest="o", help="output file")
(opt, args) = parser.parse_args()

#split to length
def split(l, n):
   return(list(splitter(l, n)))

def splitter(l, n):
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

# 1. Build sample-specific allele hash,
#    with hash[sample][chrom] pointing to a list: 
#    [pos, ref, A, B, Snptype-sample1, SNP1-sample1, SNP2-sample1, TotalCov-sample1, %A-sample1, CovA-sample1]

alleles = defaultdict(lambda: defaultdict(list))
chromlist = []

# f = open("alleles_scaleup_test.tsv", 'r')
# f = open("smol_alleles.tsv", 'r')
f = open(opt.a, 'r')

header       = f.readline().rstrip().split('\t')
snptype_cols = header[5::6]
sample_names = [x.replace("Snptype-", '') for x in snptype_cols]
print(sample_names)

linect = 0

for line in f:
    x = line.rstrip().split('\t')
    linect +=1
    if linect % 10000 == 8:
        print(linect)
    
#     # if chrom not seen before, keep track. Will iterate chroms later.
    chrom = x[0]
    if chrom not in chromlist:
        chromlist.append(chrom)
    prab = x[1:5] # Pos Ref A B = prab
    rest = split(x[5:],6) # sample-specific allele depth
    for index, name in enumerate(sample_names):
        vals = rest[index]
        if vals == ['.', '.', '.', '.', '.', '.']:
            continue
        if vals[5] == '.':
            vals[5] = 0
        if vals[4] == '.':
            vals[4] = 0
#         covs = [int(vals[3], float(vals[4])] # v2 did it this way. Want to compare.
        alleles[name][chrom].append(prab + vals)

f.close()

# 2. Go through lesion file. 
#    Build sample-specific lesion hash,
#    with hash[sample][chrom][lesion] pointing to a list: [all lesion attributes].
#    Plan to append PercALesion to this. Also append PercANotInLesion.

lesions = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
# # f = open("2019_0415_simhyb_lesions_for_testing.tsv", 'r')
# f = open("smol_indels.tsv", 'r')
f = open(opt.l, 'r')
leshead = f.readline().rstrip()

for line in f:
    x = line.rstrip().split('\t')
    sample = x[4] # updated to match real data format
    
    for names in alleles.keys():
        
        if sample in names:
            
            sample = names[:]
            chrom  = x[0]
            start  = x[1]
            end    = x[2]
            
            for z in x:
                if z == '':
                    x[x.index(x)] = '.'
                    
            lesions[sample][chrom]["{}-{}".format(start,end)] = [x[:], []]
            
f.close()

# 3. For each sampleXlesion in lesion hash,
#    ask whether loci specified in allele hash intersect the lesion. 
#    Add to sample-specific "not in lesion" hash if they don't.

notinlesion = defaultdict(list)
# 
for sample in alleles.keys():
# 
    for chrom in alleles[sample].keys():
#         
        temp = alleles[sample][chrom]
#         
        chromlesions = lesions[sample][chrom].keys()
#         print("chromlesions: {}".format(chromlesions))
#         
        for sub in temp:
            lesionhit = list(filter(lambda z: int(sub[0]) >= int(z.split('-')[0]) and int(sub[0]) <= int(z.split('-')[1]), chromlesions))
#             print("lesionhit: {}".format(lesionhit))
            if lesionhit == []:
                notinlesion[sample].append(sub)
            else:
                lesions[sample][chrom][lesionhit[0]][1].append(sub)
            if len(lesionhit) > 1:
                pass

# 4. Once all allele-specific read info is obtained for sample X lesion, determine the
#    allele ratio from the saved CallAlleles values: TotalCov and CovA.

notinvals = defaultdict(list)

# get percAnotinLes for each
for sample in notinlesion.keys():
    notLesCov  = sum(map(int, [x[7] for x in notinlesion[sample]]))
    if notLesCov == 0:
        notLesPerA = "NA"
    else:
        notLesA    = sum(map(int, [x[9] for x in notinlesion[sample]]))
        print("sample: {} notLesCov: {}".format(sample, notLesCov))
        notLesPerA = str(round((notLesA / notLesCov), 3))
    SNPNotLes = str(len(notinlesion[sample]))
    notinvals[sample] = [notLesPerA,str(notLesCov),SNPNotLes]

# o = open("test.out", 'w')
o = open(opt.o, 'w')
leshead+=('\t'.join(["\tPercALesion", "TotalCovLesion", "NumSNPLesion", "PercANotLesion", "TotalCovNotLesion","NumSNPNotLesion\n"]))
o.write(leshead)

for sample in lesions.keys():
    for chrom in lesions[sample]:
        for lesion in lesions[sample][chrom]:
            CovLesion    = sum(map(int, [x[7] for x in lesions[sample][chrom][lesion][1]]))
            if CovLesion == 0:
                percALesion = "NA"
            else:
                CovALesion   = sum(map(int, [x[9] for x in lesions[sample][chrom][lesion][1]]))
                print("sample: {} chrom: {} lesion: {} CovALesion {}".format(sample, chrom, lesion, CovALesion))
                percALesion = str(CovALesion/CovLesion)
            print(lesions[sample][chrom][lesion][1])
            SNPLes = len(lesions[sample][chrom][lesion][1])
            print("SNPles: {}".format(SNPLes))
            lesions[sample][chrom][lesion][0].extend([percALesion, str(CovLesion), str(SNPLes)])
            lesions[sample][chrom][lesion][0].extend(notinvals[sample])
            out = '\t'.join(lesions[sample][chrom][lesion][0])
            o.write(out+'\n')

o.close()