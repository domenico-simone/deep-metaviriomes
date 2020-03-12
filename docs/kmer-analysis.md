# Kmer analysis viral vs bacterial/archaeal genomes

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Set up working environment](#set-up-working-environment)
	- [NCBI genomes to be included in the analysis](#ncbi-genomes-to-be-included-in-the-analysis)
	- [Gtdb-tk installation](#gtdb-tk-installation)
- [Main analysis](#main-analysis)
	- [Calculate k-mer frequencies (jellyfish)](#calculate-k-mer-frequencies-jellyfish)
		- [Generate chunks of 5000 genomes](#generate-chunks-of-5000-genomes)
		- [Viral contigs](#viral-contigs)
		- [MAGs](#mags)
		- [NCBI genomes](#ncbi-genomes)
	- [Taxonomic assignment of MAGs from Wu et al.](#taxonomic-assignment-of-mags-from-wu-et-al)
	- [Concatenate outputs and join with taxonomic assignments](#concatenate-outputs-and-join-with-taxonomic-assignments)
	- [Taxonomic assignment with gtdb-tk](#taxonomic-assignment-with-gtdb-tk)
	- [Calculate Mean Absolute Error (MAE) to get host prediction (bacterial/archaeal vs viral)](#calculate-mean-absolute-error-mae-to-get-host-prediction-bacterialarchaeal-vs-viral)
- [Taxonomic assignment of putative host genomes not available in gtdb-tk (gtdb-tk)](#taxonomic-assignment-of-putative-host-genomes-not-available-in-gtdb-tk-gtdb-tk)
	- [Run batches](#run-batches)
	- [Summarise assignments at the relevant_rank level](#summarise-assignments-at-the-relevantrank-level)

<!-- /TOC -->

We are trying also `jellyfish`.

We also need to copy **our MAGs**, but we need first to exclude contigs which have been shown to be viral.

## Set up working environment

### NCBI genomes to be included in the analysis

The list of NCBI genomes to be included in the analysis was downloaded on March 2018 and is available in the file `data/assembly_summary.txt`. 

### Gtdb-tk installation

...

## Main analysis

### Calculate k-mer frequencies (jellyfish)

#### Generate chunks of 5000 genomes

**Try with sets of 5000 genomes.**

```bash
cd /crex2/proj/sllstore2017037/nobackup/domenico/new_analyses/aspo_viral

mkdir -p intermediate/jellyfish/NCBI_5000

python -c "
inhandle = open('data/assembly_summary.txt', 'r')
outhandle_template = 'intermediate/jellyfish/NCBI_5000/assembly_summary.links.%d'

file_n = 0
seq_n = 0
outhandle = open(outhandle_template % file_n, 'w')
for l in inhandle:
    if l.startswith('#') == False:
        l = l.split('\t')
        remoteFolder = l[19]
        fileName = l[19].split('/')[-1] + '_genomic.fna.gz'
        completePath = '/'.join([remoteFolder, fileName])
        if seq_n == 5000:
            outhandle.close()
            file_n += 1
            outhandle = open(outhandle_template % file_n, 'w')
            outhandle.write(completePath + '\n')
            seq_n = 0
        else:
            outhandle.write(completePath + '\n')
            seq_n += 1

outhandle.close()
"

ls \
intermediate/jellyfish/NCBI_5000/assembly_summary.links* > \
intermediate/jellyfish/NCBI_5000/assembly_summary.all.links
```

Run Jellyfish on batches of NCBI genome sequences.

```bash
export wdir=`pwd`
export outDir=intermediate/jellyfish/NCBI_5000

cd ${wdir}
mkdir -p ${outDir}

sbatch -p core -A snic2018-3-187 -t 50:00:00 \
-J jellyfish_NCBI_%a \
-o intermediate/jellyfish/NCBI_5000/jellyfish_NCBI_%a.out \
-e intermediate/jellyfish/NCBI_5000/jellyfish_NCBI_%a.err \
--array=1-$(wc -l < intermediate/jellyfish/NCBI_5000/assembly_summary.all.links) \
--mail-type=ALL --mail-user=domenico.simone@lnu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load jellyfish

#outDir=intermediate/jellyfish/NCBI_5000
#mkdir -p ${outDir}

genomeList=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${outDir}/assembly_summary.all.links)

echo ${genomeList}

cp ${genomeList} ${SNIC_TMP} && cd ${SNIC_TMP}

while read line; do
    echo $line
    time python -c "    
import ftputil, hashlib, sys

clean_arg = sys.argv[1].replace('ftp://ftp.ncbi.nlm.nih.gov', '')
mygenome = '/'.join(clean_arg.split('/')[1:])
genome_name = mygenome.split('/')[-1]
md5_file = '/'.join(mygenome.split('/')[:-1] + ['md5checksums.txt'])
host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')

# Download the compressed genome fasta file
# and check md5sum to ensure success

print 'Downloading genome: try 1'
host.download(mygenome, genome_name)
host.download(md5_file, 'md5checksums.txt')
md5check_exp = [l.split()[0] for l in open('md5checksums.txt', 'r').readlines() if genome_name in l][0]
md5check_fetched = hashlib.md5(open(genome_name,'rb').read()).hexdigest()
t = 1

while md5check_exp != md5check_fetched:
    host.download(mygenome, genome_name)
    md5check_fetched = hashlib.md5(open(genome_name,'rb').read()).hexdigest()
    t += 1
    print 'Downloading genome: try %d' % t

print 'Download succeeded!'
" $line
    justID=$(echo $line | awk 'BEGIN{FS="/"}{print $NF}')
    for k in 1 2 3 4 5; do
        outFile=${justID/.fna.gz/.jf${k}}
        outTab=${justID/.fna.gz/.k${k}}
        time zcat ${justID} | jellyfish count \
        -m ${k} \
        -o ${outFile} \
        -s 100M <(zcat ${justID})
        #-C ${justID}
        time jellyfish dump -c ${outFile} -o ${outTab}
    done
done < $(basename $genomeList)

cp *.jf* ${wdir}/${outDir}
cp *.k* ${wdir}/${outDir}

EOF
```

#### Viral contigs

#### MAGs

#### NCBI genomes

### Taxonomic assignment of MAGs from Wu et al.

### Concatenate outputs and join with taxonomic assignments

Many putative host genomes have no taxonomy. We should get their AC, download them and perform taxonomic assignment.

### Taxonomic assignment with gtdb-tk

### Calculate Mean Absolute Error (MAE) to get host prediction (bacterial/archaeal vs viral)

## Taxonomic assignment of putative host genomes not available in gtdb-tk (gtdb-tk)

### Run batches

### Summarise assignments at the relevant_rank level

```python
a = open('results/viral_contigs_host_prediction/MAE/viral_contigs_host_assignment_k4.181114.tsv', 'r')

# dict of viral contigs : relevant ranks
D = {}

for l in a:
    l = l.split()
    if len(l) > 3:
        if l[0] not in D:
            D[l[0]] = [set([l[-1]]), set()]
        else:
            D[l[0]][0].add(l[-1])
        if l[2][0] == "P":
            D[l[0]][1].add(l[-1])

b = open('results/viral_contigs_host_prediction/MAE/viral_contigs_host_assignment_k4.181114.summary.tsv', 'w')
for k in D:
    summary = ','.join(list(D[k][0]))
    summary_MAG = ','.join(list(D[k][1]))
    b.write('\t'.join([k, summary_MAG, summary]) + '\n')

b.close()
```