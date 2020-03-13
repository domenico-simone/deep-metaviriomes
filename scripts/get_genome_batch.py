import os, glob, sys

def get_genomeID(wholepath):
    return wholepath.split('/')[-1].split('.')[0]

# create output folder
try:
    os.mkdir("input_batches")
except:
    pass

# list of genome files
L = glob.glob("data/genomes/aspo_MAGs/*.fa")

# counter of batch files
b = 1
# counter of genome files
g = 1
# batch_file
bfile = open('input_batches/genomes_batch.%d' % b, 'w')

for f in L:
    if g%40 == 0:
        b += 1
        bfile = open('input_batches/genomes_batch.%d' % b, 'w')
    #    bfile.write('\t'.join([f, get_genomeID(f)])+'\n')
    #else:
    bfile.write('\t'.join([f, get_genomeID(f)])+'\n')
    g += 1

bfile.close()
