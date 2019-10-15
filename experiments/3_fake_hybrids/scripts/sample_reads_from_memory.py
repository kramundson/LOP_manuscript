#!/share/comailab/kramundson/miniconda3/envs/ximena/bin/python

import argparse
import gzip
import random
import sys
import time

"""
Using raw fastq reads, generate pseudorandom subsets of single-end reads from 2 parents.
Does not work on paired-end reads!

Use this for simulating hybrids in cases where I don't have them.

Script takes in gzip fastq and writes uncompressed fastq.
A later script compresses each read set in parallel for faster execution.

USAGE: python script.py parent1_fastq parent2_fastq seed_start number_of_samples parent1_readnum parent2_readnum

Example command, which subsamples from files SRR6123041.fastq.gz SRR6123042.fastq.gz
to make 2 subsets composed of 100 reads from each file:
python sample_reads_from_memory.py -p1 SRR6123041.fastq.gz -p2 SRR6123042.fastq.gz -s 1 -n 2 -p1n 100 -p2n 100

In the first subset, 100 reads are pseudo-randomly subsampled from SRR6123041.fastq.gz
using a seed of 1, and 100 reads are subsampled from SRR6123042.fastq.gz using a seed of 3.
Read subsamples from each file are then concatenated and written to disk.

In the second subset, 100 reads are subsampled from SRR6123041.fastq.gz using a seed of 2
and 100 reads are subsampled from SRR6123042.fastq.gz using a seed of 4. Read subsamples
from each file are then are then concatenated and written to disk.

Uniparental test case, in which 200 reads drawn from SRR6123041.fastq.gz only:
python sample_reads_from_memory.py -p1 SRR6123041.fastq.gz -s 1 -n 2 -p1n 200


For the real deal, makes 1000 subsamples of 3023119 reads from high-coverage libraries of
LOP868 and PL4 in a 1:2 ratio:
python sample_reads_from_memory.py -p1 SRR6123032_1.fastq.gz -p1n 2015413 -p2 SRR6123031_1.fastq.gz -p2n 1007706 -s 1 -n 1000

Run failed at subsampling 851 for reasons unknown. Finish off by deleting sample 851
and restarting
rm 851_SRR6123032_SRR6123031_1851.fastq
python sample_reads_from_memory.py -p1 SRR6123032_1.fastq.gz -p1n 2015413 -p2 SRR6123031_1.fastq.gz -p2n 1007706 -s 851 -n 149


Second real deal: makes 1000 subsamples of 3023119 from high coverage libraries of LOP868
and IVP101 in a 1:2 ratio:
python sample_reads_from_memory.py -p1 SRR6123032_1.fastq.gz -p1n 2015413 -p2 SRR6123183_1.fastq.gz -p2n 1007706 -s 2001 -n 1000


Third real deal: take 1000 subsamples of 3023119 reads only from LOP868:
python sample_reads_from_memory.py -p1 SRR6123032_1.fastq.gz -p1n 3023119 -s 4001 -n 1000
"""

def gzip_reads_to_list(readfile):
    start = time.time()
    reads = []
    with gzip.open(readfile, 'rt') as f:
        while 1:
            newread = f.readline()
            if newread == "":
                break
            newread += f.readline() # seq
            newread += f.readline() # plus
            newread += f.readline() # qual
            
            reads.append(newread)
            
    finish = time.time()
    
    print("Read in took {} seconds".format(finish-start))
    print("Size of list: {} bytes".format(sys.getsizeof(reads)))
    print("List contains {} reads".format(len(reads)))
    
    return reads
    
parser = argparse.ArgumentParser(description="subsample reads from fastq")
parser.add_argument("-p1", "--parent1", dest="p1file", default=False, help="Parent1 fastq reads, gzipped")
parser.add_argument("-p1n", "--parent1n", dest="p1n", type=int, help="Number of p1 reads per draw")
parser.add_argument("-p2", "--parent2", dest="p2file", default=False, help="Parent2 fastq reads, gzipped")
parser.add_argument("-p2n", "--parent2n", dest="p2n", type=int, help="Number of p2 reads per draw")
parser.add_argument("-s", "--seed", dest="start_seed", type=int, help="starting seed")
parser.add_argument("-n", "--ndraws", dest="ndraws", type=int, help="Number of draws")

args = parser.parse_args()

start = time.time()

if not args.p1file:
    sys.exit("Provide at least one fastq file on the command line")

p1_reads = gzip_reads_to_list(args.p1file)

if args.p2file:
    p2_reads = gzip_reads_to_list(args.p2file)

print("Done reading in parents, now subsampling and writing out")

# subsample reads from p1_reads
# uses sequential seeds
# be sure that seeds are not being reused! may want to code this in!
for seed in range(args.start_seed, args.start_seed+args.ndraws+1):
    # write out to disk
    p1 = args.p1file.split('.')[0].split("/")[-1].split("_")[0]
    
    # set seed for sampling p1
    random.seed(seed)
    print("Sampling {} reads from parental reads {} using seed {}".format(args.p1n, p1, seed))
    
    # sample reads from p1
    out=random.sample(p1_reads, args.p1n) # random.sample() returns a list
    
    # set seed for sampling p2
    if args.p2file:
        p2 = args.p2file.split('.')[0].split("/")[-1].split("_")[0]
        p2_seed = seed + args.ndraws
        random.seed(p2_seed)
        print("Sampling {} reads from parental reads {} using seed {}".format(args.p2n, p2, p2_seed))
        # sample reads from p2
        out += random.sample(p2_reads, args.p2n)
        ofh = "data/reads/{}_{}_{}_{}.fastq".format(seed, p1, p2, p2_seed)
    else:
        ofh = "data/reads/{}_{}_uniparental.fastq".format(seed,p1)
        
    o = open(ofh, 'w')
    o.write("".join(out))
    o.close()

finish = time.time()
print("Script took {} seconds".format(finish-start))