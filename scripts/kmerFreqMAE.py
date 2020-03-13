#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys, os, glob, time
from sys import getsizeof
from collections import OrderedDict

def readJellyOut(infile):
    b = pd.read_table(infile, header=None, sep = " ")
    b.columns = ['kmer', 'count']
    # b = open(infile, 'r').readlines()
    # b = [[l.split()[0], int(l.split()[1])] for l in b]
    #b = np.array(b)
    return b

def calcFreqMAE(x, y):
    ## https://raw.githubusercontent.com/derekjhyang/python-mae/master/mae.py
    #return sum(map(lambda t:abs(float(t[0]/sum(x)-t[1]/sum(y))),zip(x, y)))/len(x)
    return sum(map(lambda t:abs(float(t[0])/sum(x)-float(t[1])/sum(y)),zip(x, y)))/len(x)

k = int(sys.argv[1])
batch_file = sys.argv[2]

# default:
output_basedir = "intermediate/jellyfish"
# test
# output_basedir = "intermediate/jellyfish_test"

try: os.makedirs("%s/viral_mae/k%d" % (output_basedir, k))
except: pass

try: os.makedirs("%s/viral_mae_OrderedDict/k%d" % (output_basedir, k))
except: pass

# version that re-reads all profiles and joins kmer tables first to ensure
# we're comparing the same kmer counts
#for k in kmers:
outhandle = open("%s/viral_mae/k%d/kmer_mae.k%d.%s.out" % (output_basedir, k, k, batch_file.split("/")[-1]), "w")
print "k: %d" % k

microbialSpectrum_files = glob.glob("intermediate/jellyfish/NCBI_5000/*.k%d" % k) + glob.glob("intermediate/jellyfish/MAGs/*.k%d" % k)

viralSpectrum_files = [i.strip() + ".k%d" % k for i in open(batch_file, 'r').readlines()]
for v in viralSpectrum_files[:5]:
    print "viral contig: %s" % v
    D_coeff = OrderedDict()
    viralSpectrum = readJellyOut(v)
    for x,m in enumerate(microbialSpectrum_files[:5]):
        print m
        if x%10000 == 0:
            print x
        microbialSpectrum = readJellyOut(m)
        joinedSpectra = pd.merge(viralSpectrum, microbialSpectrum, on = "kmer", how = "outer").fillna(0)
        #print joinedSpectra
        coeff = calcFreqMAE(joinedSpectra['count_x'], joinedSpectra['count_y'])
        print coeff
        ID = m.split('/')[-1].split('.')[0]
        D_coeff[ID] = coeff
    
    D_coeff = OrderedDict(sorted(D_coeff.items(), key = lambda t: t[1], reverse = True))
    #print D_coeff
    # The file here produced can be open as in https://stackoverflow.com/a/34798344/10554322
    dict_outfile = open("%s/viral_mae_OrderedDict/k%d/%s_orderedDict.out" % (output_basedir, k, v.split("/")[-1]), 'w')
    dict_outfile.write(str(D_coeff))
    dict_outfile.close()
    # Write out results with MAE <= 10^-3
    best_res = filter(lambda x: x[1] < 0.001, D_coeff.items())
    for i in best_res:
        outhandle.write('\t'.join([v.split("/")[-1], str(i[1]), str(i[0])]) + '\n')

outhandle.close()
