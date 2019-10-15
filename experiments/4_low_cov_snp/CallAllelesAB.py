# Isabelle Henry, last modified June 2012
# UC Davis Genome Center

# This script takes a parsed mpileup file as input (samtools) that contains allelic calls for various individuals.
# The mpileup comes from aligning each of these samples to the reference
# genome.

# another file is provided that contains the SNP positions between the two parental alleles and the alleles themselves

# This script compares the alleles in the two parental genotypes and assign an allele for the samples
# The criteria are as follows:
# - if only one allele type is found, and that allele corresponds to one of the two parental allele, that allele is called.
# - if there are two allele types and the least frequent allele represents < 10% of the calls, it is treated as if only the most common type was present
# - if there are two allele types and the least frequent allele represents > 10% of the allele calls, and the two alleles
# correspond to the two parents, the allele is called Het.
# - in all other cases, the genotype is called "na"

# USAGE: CallAllelesDN.py parsed_mpileup_marker.txt output.txt SNPfile.txt

# The script outputs the Chromosome, position and ref base, followed by the A and B alleles (from the SNP file).
# Next, the script outputs 6 columns per sample:
# SNP type (1 = homozygous and 2 = het)
# SNP1: Most common allele
# SNP2: Second allele
# Total cov: Total coverage at that position
# %A: 1 if homozygous A, 0 if homozygous B, 0.5 if heterozygous
# CovA: coverage of the A allele

import sys

mpup = open(sys.argv[1]) #mpileup file
o = open(sys.argv[2],'w') #output file

## make a list of SNP position and their base call in the two parental genotypes

SNPfile = open(sys.argv[3]) #snp file
SNPs = {}

SNPfile.readline()

for mine in SNPfile:
    m = mine.split('\t')
    SNPChrom = m[0]
    SNPPos = m[1]
    SNPRef = m[2].upper()
    SNP_A = m[3].upper()
    SNP_B = m[4].split('\n')[0].upper()
    SNPMega = int(SNPPos)/1000000
    if SNPChrom not in SNPs:
        SNPs[SNPChrom] = {}
    if SNPMega not in SNPs[SNPChrom]:
        SNPs[SNPChrom][SNPMega] = {}
    SNPs[SNPChrom][SNPMega][SNPPos]=[SNPRef,SNP_A,SNP_B]
SNPfile.close()

## Going through the mpileup now

indices = []
totalcount = 0
odd1 = 0
odd2 = 0

line = mpup.readline()
l = line.split()
if line[0] == 'C':
    count = 0
    header = 'Chrom\tPos\tRef\tA\tB'
    for m in l:
        if m.split('-')[0] == 'Cov':
            indices.append(count)
            lib = m[4:]
            header += '\tSnptype-'+lib+'\tSNP1-'+lib+'\tSNP2-'+lib+'\tTotalCov-'+lib+'\t%A-'+lib+'\tCovA-'+lib
        count += 1
    o.write(header+'\n')
    print indices

for line in mpup:
    l = line.split()
    chrom = l[0]
    pos = l[1]
    mega = int(pos)/1000000
    ref = l[2].upper()
    text = ''
    if chrom in SNPs:
        if mega in SNPs[chrom]:
            if pos in SNPs[chrom][mega]:
                a = SNPs[chrom][mega][pos][1].upper()
                b = SNPs[chrom][mega][pos][2].upper()
                refb = SNPs[chrom][mega][pos][0].upper()
                if refb != ref:
                    print "ref problem"
                    break
                else:
                    text = chrom+'\t'+str(pos)+'\t'+ref+'\t'+a+'\t'+b

                for i in indices:
                    Cov = l[i]
                    if Cov == '.':
                        text += 6 * '\t.'
                    else:
                        totalcount += 1
                        snp1 = l[i-3]
                        snp2 = l[i-2]
                        snp3 = l[i-1]
                        if snp3 != '.':
                            Snptype = 3
                        elif snp2 != '.':
                            Snptype = 2
                        else:
                            Snptype = 1
                        Snp1 = snp1.split('_')[0]
                        # print(l) # debugging
                        # print(Snp1)
                        # print(i)
                        Cov1 = int(round(float(snp1.split("_")[1])*float(Cov)/100)) # is throwing index error
                        if snp2 == '.':
                            Snp2 = '.'
                        elif float(snp2.split('_')[1]) <= 10:
                            Snp2 = '.'
                        else:
                            Snp2 = snp2.split('_')[0]
                            Cov2 = int(round(float(snp2.split('_')[1])*float(Cov)/100))
                        if (Snp1 == a and Snp2 == b) or (Snp1 == b and Snp2 == a):
                            if Snp1 == a:
                                All = 0.5
                                CovA = Cov1
                            elif Snp1 == b:
                                All = 0.5
                                CovA = Cov2
                        elif Snp2 != '.':
                            All = '.'
                            odd2 += 1
                            CovA = '.'
                        elif Snp1 == b:
                            All = 0
                            CovA = 0
                        elif Snp1 == a:
                            All = 1
                            CovA = Cov
                        else:
                            All = '.'
                            odd1 += 1
                            CovA = '.'
                        text += '\t'+str(Snptype)+'\t'+str(Snp1)+'\t'+str(Snp2)+'\t'+str(Cov)+'\t'+str(All)+'\t'+str(CovA)
                text += '\n'
                o.write(text)
    
o.close()
mpup.close()

print totalcount
print odd1
print odd2
