# Isabelle Henry, last modified June 2012
# UC Davis Genome Center

# This script takes a file of genotypes (output of CalLAllelesAB script) and outputs a mean per chromosome or per bin for each sample

# USAGE: bin-by-genotype.py alleles output.txt binsize snpfile

# where alleles is the output of the CallAllelesAB script
# output is the desired name of the output file
# binsize is the desired binsize in bps

import sys

alleles = open(sys.argv[1],'r')
o = open(sys.argv[2],'w')
binna = int(sys.argv[3])

# Block added KRA 12/22/16
SNP = open(sys.argv[4],'r') # add to record number of parent-informative SNPs per bin
maxBin = {} # add to record number of parent-informative SNPs per bin

# Going through SNP file, recording number of parent-informative SNPs per bin
SNP.readline()
for line in SNP:
    l = line.split()
    chrom = l[0]
    SNPpos = int(l[1])/binna
    if len(str(SNPpos)) < 7:
        SNPbinS = '0'*(7-len(str(SNPpos)))+str(SNPpos) # add required number of 0's to end
    else:
        SNPbinS = str(SNPpos) # if no 0's needed, then just add str.
    try:
        BinNum = int(chrom.split('chr')[1]+str(SNPbinS)) #tack on chromosome number
    except:
        continue # note, skips ChrUn
    if BinNum not in maxBin:
        maxBin[BinNum] = 1
    else:
        maxBin[BinNum] += 1
SNP.close()
# end KRA added block




## Going through the allele file

totalcount = 0

line = alleles.readline()
l = line.split()

Libs = {}
Bins = {}

if line[0] == 'C':
    count = 0
    l = line.split('\t')
    for m in l:
        if m.split('-')[0] == 'TotalCov':
            lib = str(m[7:])
            Libs[lib,count] = {}
        count += 1
    print Libs

for line in alleles:
    l = line.split()
    #if l[0][0:2] == "Ch": #case-sensitive, comment out
    if l[0][0:2] == "ch":
        chrom = l[0]
#         print(chrom)
        pos = int(l[1])
        bins = int(pos/binna) # modified in fruitless attempt at debugging
        if len(str(bins)) < 7:
            binS = '0'*(7-len(str(bins)))+str(bins)
        else:
            binS = str(bins)
#             print(binS)
        
        # KRA added try:except block 12/22/16 to accommodate case sensitivity. Note: ChrUn not analyzed here.
        try:
#             print "try runs" # try running 4 times
            BinNum = int(chrom.split("chr")[1]+str(binS))
        except:
#             print "except runs" #only runs four times?
            continue
        # end KRA added block
        
        if BinNum not in Bins:
            Bins[BinNum] = [chrom,bins*binna,(bins*binna)+binna,maxBin[BinNum]] #KRA modified line to count parent-informative SNPs per bin
#             print("Bins[BinNum]: "+ str(Bins[BinNum])) # debugging
            for Lib in Libs:
                Libs[Lib][BinNum] = {}
                Libs[Lib][BinNum]['Cov'] = 0
                Libs[Lib][BinNum]['ExpACov'] = 0
                Libs[Lib][BinNum]['ObsACov'] = 0
        for i in Libs:
            if l[i[1]+2] != '.':
                # add the cov
                Libs[i][BinNum]['Cov'] += int(l[i[1]])
                # add the expected number of instance of the A allele (genotype x cov)
                Libs[i][BinNum]['ExpACov'] += int(float(l[i[1]])*float(l[i[1]+1]))
                # add the observed number of instances of the A allele
                Libs[i][BinNum]['ObsACov'] += int(l[i[1]+2])

alleles.close()
# print(Bins) # debugging., works fine
#print(maxBin) # debugging, works fine

sortedbins = Bins.keys()
sortedbins.sort()

sortedlibs = Libs.keys()
sortedlibs.sort()

header = "Chrom\tStart\tEnd\tMax"

for Sample in sortedlibs:
    header += "\t"+Sample[0]+"-Cov"

for Sample in sortedlibs:
    header += "\t"+Sample[0]+"-Obs%A"

for Sample in sortedlibs:
    header += "\t"+Sample[0]+"-Calc%A"

o.write(header+'\n')

for i in sortedbins:
    Chrom = str(Bins[i][0])
    Start = str(Bins[i][1])
    End = str(Bins[i][2])
    Max = str(Bins[i][3])
    text = Chrom + '\t'+ Start + '\t' + End + '\t' + Max
    
    for Sample in sortedlibs:
        coverage = Libs[Sample][i]['Cov']
        text += '\t'+str(coverage)
    
    for Sample in sortedlibs:
        coverage = Libs[Sample][i]['Cov']
        if coverage > 0:
            ObsPerA = 100*(Libs[Sample][i]['ObsACov'])/(Libs[Sample][i]['Cov'])
        else:
            ObsPerA =  "."
        text += '\t'+str(ObsPerA)

    for Sample in sortedlibs:
        coverage = Libs[Sample][i]['Cov']
        if coverage > 0:
            CalcPerA = 100*(Libs[Sample][i]['ExpACov'])/(Libs[Sample][i]['Cov'])
        else:
            CalcPerA =  "."
        text += '\t'+str(CalcPerA)

    o.write(text+'\n')

o.close()