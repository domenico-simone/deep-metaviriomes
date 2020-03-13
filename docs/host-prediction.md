# Host prediction for viral contigs

<!-- TOC START min:2 max:5 link:true asterisk:false update:true -->
- [Set up working environment](#set-up-working-environment)
    - [Gtdb-tk installation](#gtdb-tk-installation)
- [Data included in this protocol](#data-included-in-this-protocol)
    - [Viral contigs](#viral-contigs)
    - [NCBI genomes](#ncbi-genomes)
    - [Äspö MAGs not in NCBI](#äspö-mags-not-in-ncbi)
- [Main analysis](#main-analysis)
    - [Calculate k-mer frequencies (jellyfish)](#calculate-k-mer-frequencies-jellyfish)
        - [NCBI bacterial/archaeal genomes](#ncbi-bacterialarchaeal-genomes)
        - [Viral contigs](#viral-contigs-1)
        - [Äspö MAGs](#äspö-mags)
    - [what](#what)
    - [Taxonomic assignment of NCBI genomes and Äspö MAGs](#taxonomic-assignment-of-ncbi-genomes-and-äspö-mags)
        - [Taxonomic assignment of NCBI genomes missing from GTDB](#taxonomic-assignment-of-ncbi-genomes-missing-from-gtdb)
        - [Taxonomic assignment of MAGs from Wu et al.](#taxonomic-assignment-of-mags-from-wu-et-al)
    - [Concatenate outputs and join with taxonomic assignments](#concatenate-outputs-and-join-with-taxonomic-assignments)
    - [Taxonomic assignment with gtdb-tk](#taxonomic-assignment-with-gtdb-tk)
    - [Calculate Mean Absolute Error (MAE) to get host prediction (bacterial/archaeal vs viral)](#calculate-mean-absolute-error-mae-to-get-host-prediction-bacterialarchaeal-vs-viral)
- [Taxonomic assignment of putative host genomes not available in gtdb-tk (gtdb-tk)](#taxonomic-assignment-of-putative-host-genomes-not-available-in-gtdb-tk-gtdb-tk)
    - [Run batches](#run-batches)
    - [Re-merge MAE results](#re-merge-mae-results)
    - [Summarise assignments at the relevant_rank level](#summarise-assignments-at-the-relevant_rank-level)
<!-- TOC END -->


We are trying also `jellyfish`.

We also need to copy **our MAGs**, but we need first to exclude contigs which have been shown to be viral.

## Set up working environment

This analysis has been run on a computing cluster with a SLURM job scheduler. Access to resources managed by SLURM usually requires an active account where computing time is charged. In order to replicate the analyses described in this repo, you need to set the following variable which identifies the account on the cluster.

```bash
export SLURMaccount="snic2018-3-187"
```

### Gtdb-tk installation

Use **gtdb-tk**. Install it as conda environment, as described here https://bitbucket.org/scilifelab-lts/m_dopson_1701/wiki/Installation.

## Data included in this protocol

### Viral contigs

The dataset included here is available at NCBI with accession numbers:MT141153-MT145203.

### NCBI genomes

The list of NCBI genomes to be included in the analysis was downloaded on March 2018 and is available in the file `data/assembly_summary.txt`. 

### Äspö MAGs not in NCBI

Metagenome assembled genomes (hereafter named **MAGs**) binned from ...

## Main analysis

### Calculate k-mer frequencies (jellyfish)

#### NCBI bacterial/archaeal genomes

Generate chunks of 5000 genomes (`intermediate/jellyfish/NCBI_5000/assembly_summary.links.*`) and save a list of them in the file `intermediate/jellyfish/NCBI_5000/assembly_summary.all.links`.

```bash
#cd /crex2/proj/sllstore2017037/nobackup/domenico/new_analyses/aspo_viral

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

Run Jellyfish on batches of NCBI bacterial/archaeal genome sequences

```bash
export wdir=`pwd`
export outDir=intermediate/jellyfish/NCBI_5000

cd ${wdir}
mkdir -p ${outDir}

sbatch -p core -A ${SLURMaccount} -t 50:00:00 \
-J jellyfish_NCBI_%a \
-o intermediate/jellyfish/NCBI_5000/jellyfish_NCBI_%a.out \
-e intermediate/jellyfish/NCBI_5000/jellyfish_NCBI_%a.err \
--array=1-$(wc -l < intermediate/jellyfish/NCBI_5000/assembly_summary.all.links)<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load jellyfish

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
    for k in 4; do
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

Get separate lists of viral contigs to have 4 separate jobs, in order to distribute the computational burden.
Create output directory `intermediate/jellyfish/viral`.

```bash
ls data/genome_viral/P911_10*.fa | awk 'NR<1013' > data/genome_viral/contigList.1
ls data/genome_viral/P911_10*.fa | awk 'NR>1012&&NR<2026' > data/genome_viral/contigList.2
ls data/genome_viral/P911_10*.fa | awk 'NR>2025&&NR<3039' > data/genome_viral/contigList.3
ls data/genome_viral/P911_10*.fa | awk 'NR>3038' > data/genome_viral/contigList.4

mkdir -p intermediate/jellyfish/viral
```

Run one job for each contigList file

```bash
for f in $(ls ./data/genome_viral/contigList.*); do
    export f=$f
    serial=$(basename $f | awk 'BEGIN{FS="."}{print $NF}')
    #echo $f, $q
    sbatch -p core -n 10 -A ${SLURMaccount} -t 4:00:00 \
    -J jellyFish.${serial} -o logs/jellyFish.${serial}.out -e logs/jellyFish.${serial}.err \
    --mail-type=ALL --mail-user=domenico.simone@lnu.se<<'EOF'
#!/bin/bash

module load jellyfish
outDir=intermediate/jellyfish/viral

for k in 4; do
    for contigFile in $(cat ${f}); do
        b=$(basename ${contigFile})
        outFile=$(echo ${b/.fa/.jf${k}})
        outTab=$(echo ${b/.fa/.k${k}})
        jellyfish count \
        -m${k} \
        -o ${outDir}/${outFile} \
        -s 100M \
        -t 10 \
        ${contigFile}
        jellyfish dump -c ${outDir}/${outFile} -o ${outDir}/${outTab}
    done
done

EOF
done
```

#### Äspö MAGs

```bash
ls /crex2/proj/sllstore2017037/nobackup/shared_data/aspo/metagenomes/scilife_metagenomes/data_papers_planktonic_biofilm/bins50/genomes/*/*/*.fa > data/MAG_list

# create output directory for jellyfish
mkdir -p intermediate/jellyfish/MAGs

export f=data/MAG_list
#serial=$(basename $f | awk 'BEGIN{FS="."}{print $NF}')
#echo $f, $q
sbatch -p core -n 10 -A ${SLURMaccount} -t 4:00:00 \
-J jellyFish.MAG -o jellyFish.MAG.out -e jellyFish.MAG.err \
--mail-type=ALL --mail-user=domenico.simone@lnu.se<<'EOF'
#!/bin/bash

module load jellyfish
outDir=intermediate/jellyfish/MAGs

for k in 4; do
    for contigFile in $(cat ${f}); do
        b=$(basename ${contigFile})
        outFile=$(echo ${b/.fa/.jf${k}})
        outTab=$(echo ${b/.fa/.k${k}})
        jellyfish count \
        -m ${k} \
        -o ${outDir}/${outFile} \
        -s 100M \
        -t 10 \
        ${contigFile}
        jellyfish dump -c ${outDir}/${outFile} -o ${outDir}/${outTab}
    done
done

EOF
```

#### Create list of batches of jellyfish results for viral contigs 

Create list of batches

```python
# python

import glob, os

A = glob.glob("intermediate/jellyfish/viral/*.k4")
A = [i.replace(".k1", "") for i in A]

outdir = "intermediate/jellyfish/viral/batches"
try: os.mkdir(outdir)
except: pass

c = 1
t = 1
for i in A:
    outhandle = open('%s/viral_batch.%s' % (outdir, str(t).zfill(4)), 'w')
    if c%5 == 0:
        t += 1
        outhandle.close()
        outhandle = open('%s/viral_batch.%s' % (outdir, str(t).zfill(4)), 'w')
    outhandle.write(i + '\n')
    c += 1
```

```bash
ls intermediate/jellyfish/viral/batches/viral_batch.* > intermediate/jellyfish/viral/batches/viral_batch_list
```


## Calculate MAE of kmer frequencies to get host prediction

As in doi:10.1038/nature19366.

The script `kmerFreqMAE.py` (in the `scripts` folder) performs calculation of mean absolute error between all viral contig kmer frequency profiles and NCBI/Aspo genome kmer frequency profiles.

Raw results are stored in folder `mkdir -p intermediate/jellyfish/viral_mae_OrderedDict/k4`.

Best results (MAE < 0.001) are stored in folder `mkdir -p intermediate/jellyfish/viral_mae/k4`.

```bash
# Create output directories

mkdir -p intermediate/jellyfish/viral_mae_OrderedDict/k4
mkdir -p intermediate/jellyfish/viral_mae/k4
mkdir -p logs/viral_mae

export k=4

sbatch -n1 -p core -A snic2018-3-187 -t 2-00:00:00 \
--array=1-$(wc -l < intermediate/jellyfish/viral/batches/viral_batch_list) \
-J kmer_freq_mae.k${k}_%a -e logs/viral_mae/kmer_freq_mae.k${k}_%a.err -o logs/viral_mae/kmer_freq_mae.k${k}_%a.out \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

python scripts/kmerFreqMAE.py ${k} $(sed -n "$SLURM_ARRAY_TASK_ID"p intermediate/jellyfish/viral/batches/viral_batch_list)

EOF
```

### what

```python
import pandas as pd
import glob, os

def clean(NCBI_AC):
    a = NCBI_AC
    try:
        NCBI_AC = NCBI_AC.split('_')[1:]
        NCBI_AC[-1] = NCBI_AC[-1].split('.')[0]
        NCBI_AC = '_'.join(NCBI_AC)
    except:
        NCBI_AC = a
        #print "Cleaning yielded an error. Got %s, AC is %s" % (a, NCBI_AC)
    return NCBI_AC

def trick(cell):
    if "GCA" in cell:
        cell = cell.replace("GCA", "GCF")
    elif "GCF" in cell:
        cell = cell.replace("GCF", "GCA")
    else:
        cell = cell
    return cell

def getRelevantRank(taxonomy):
    relevant_rank = ""
    #taxonomy = taxonomy.split(";")
    try:
        if "Proteobacteria" in taxonomy:
            taxonomy = taxonomy.split(";")
            for i in taxonomy:
                if i.startswith("c__"):
                    relevant_rank = i[3:]
        else:
            taxonomy = taxonomy.split(";")
            for i in taxonomy:
                if i.startswith("p__"):
                    relevant_rank = i[3:]
    except:
        print "Taxonomy %s yielded an error" % taxonomy
    if relevant_rank == "":
        print "Taxonomy %s returned an empty relevant rank" % taxonomy
    return relevant_rank

def taxonomic_assignment(kmer_value):
    kmer_value = int(kmer_value)
    path = r'intermediate/jellyfish/viral_cc/k%d' % kmer_value
    allFiles = glob.glob(path + "/*viral_batch*.out")
    frame = pd.DataFrame()
    list_ = []
    for file_ in allFiles:
        df = pd.read_table(file_, header=None, sep = "\t")
        list_.append(df)
    
    frame = pd.concat(list_)
    frame.columns = ['viral_contig.kmer', 'top_hits_cc', 'top_hits']
    
    # read NCBI taxonomic assignments
    NCBI_taxonomy = pd.read_table("data/bac_taxonomy_r86.tsv")
    NCBI_taxonomy.columns = ["AC", "taxonomy"]
    NCBI_taxonomy['AC'] = NCBI_taxonomy['AC'].apply(clean)
    
    # since some AC are shifted from GCA to GCF,
    # create a fake version of the NCBI_taxonomy table
    
    NCBI_taxonomy_fake = NCBI_taxonomy.copy()
    NCBI_taxonomy_fake['AC'] = NCBI_taxonomy_fake['AC'].apply(lambda x: trick(x))
    
    # read MAG taxonomic assignments
    MAG_taxonomy_files = glob.glob("../gtdb-tk_test/gtdbtk_output_batch.*/*classification_pplacer.tsv")
    #print MAG_taxonomy_files
    mframe = pd.DataFrame()
    mlist_ = []
    for file_ in MAG_taxonomy_files:
        mdf = pd.read_table(file_, header=None, sep = "\t")
        mlist_.append(mdf)
    
    mframe = pd.concat(mlist_)
    mframe.columns = ["AC", "taxonomy"]
    
    # rbind taxonomy tables
    all_taxonomy = NCBI_taxonomy.append(mframe).append(NCBI_taxonomy_fake)
    
    # join!
    frame_taxonomy = pd.merge(frame, all_taxonomy, left_on = "top_hits", right_on = "AC", how = "left")
    
    # get phylum or class, in case of Proteobacteria
    frame_taxonomy['taxonomy_short'] = frame_taxonomy['taxonomy'].apply(lambda x: getRelevantRank(x))
    frame_taxonomy.to_csv(path_or_buf = "results/viral_contigs_host_assignment_k%d.tsv" % kmer_value, sep = "\t", index = False)
    #return frame_taxonomy

taxonomic_assignment(4)
```

### Taxonomic assignment of predicted hosts

For most NCBI genomes used in this analysis taxonomic assignments are already available through GTDB (file `data/bac_taxonomy_r86.tsv`). For the ones with no taxonomic assignment available through GTDB, we'll compute it ourselves with gtdb-tk. We need also to assign taxonomy to the Äspö MAGs.

### Taxonomic assignment of NCBI genomes and Äspö MAGs

#### Taxonomic assignment of NCBI genomes missing from GTDB

Get list of NCBI genomes with missing taxonomic assignment. 

```python
# python

import pandas as pd
import glob, os, sys

unknown_genomes = set()

for outf in glob.glob('results/viral_contigs_host_assignment_k*.tsv'):
    a = open(outf, 'r')
    for i in a:
        i = i.split("\t")
        if i[3] == "":
            unknown_genomes.add(i[2])

#try:
os.makedirs("intermediate/unknown_genomes")
ppp = open("intermediate/unknown_genomes/unknown_genomes_list", 'w')
for i in unknown_genomes:
    ppp.write(i+'\n')

ppp.close()
os.system("grep -hf intermediate/unknown_genomes/unknown_genomes_list data/assembly_summary.links.* > intermediate/unknown_genomes/unknown_genomes_list.links")
#except: sys.exit("Couldn't create dir.")
#len(unknown_genomes)
```

#### Taxonomic assignment of MAGs from Wu et al.

### Concatenate outputs and join with taxonomic assignments

Many putative host genomes have no taxonomy. We should get their AC, download them and perform taxonomic assignment.

### Taxonomic assignment with gtdb-tk

Taxonomic assignment with **gtdb-tk**

```bash
# activate conda env if not active yet
conda activate /crex/proj/sllstore2017037/nobackup/domenico/gtdb-tk

# generate log folder
mkdir -p logs/gtdbtk/unknown_genomes_2

# generate genome batches
ls intermediate/unknown_genomes_2/*.fna.gz > intermediate/unknown_genomes_2/unknown_genomes.ls
cd intermediate/unknown_genomes_2 && split -l 200 --numeric-suffixes unknown_genomes.ls unknown_genomes.ls.batch && cd -

for i in $(ls intermediate/unknown_genomes_2/unknown_genomes.ls.batch*); do
    export batchPath=${i}
    export batchFile=$(basename ${i})
    export export_dir=`pwd`
sbatch -A snic2017-7-182 -p node -t 30:00:00 \
-J gtdbtk_unknown.${batchFile} -o logs/gtdbtk/unknown_genomes_2/gtdbtk_unknown.${batchFile}.out -e logs/gtdbtk/unknown_genomes_2/gtdbtk_unknown.${batchFile}.err \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

for j in $(cat ${batchPath}); do
    cp ${j} ${SNIC_TMP}
done

cd ${SNIC_TMP}

for i in $(ls *fna.gz); do
    zcat ${i} > ${i/.fna.gz/.fna}
done

ls *.fna | awk 'BEGIN{FS=".";OFS="\t"}{print $0, $1}' > unknown_genomes.batchfile

gtdbtk classify_wf --cpus 20 --min_perc_aa 50 --debug --batchfile unknown_genomes.batchfile --out_dir unknown_genomes.${batchFile}_gtdbtk_out -x fna

rsync -av unknown_genomes.${batchFile}_gtdbtk_out ${export_dir}/intermediate/unknown_genomes_2

EOF
done
```

### Calculate Mean Absolute Error (MAE) to get host prediction (bacterial/archaeal vs viral)

## Taxonomic assignment of putative host genomes not available in gtdb-tk (gtdb-tk)

### Run batches


### Re-merge MAE results

Re-merge MAE analysis with taxonomic assignments, adding the new ones

```python
import pandas as pd
import glob, os

def clean(NCBI_AC):
    a = NCBI_AC
    try:
        NCBI_AC = NCBI_AC.split('_')[1:]
        NCBI_AC[-1] = NCBI_AC[-1].split('.')[0]
        NCBI_AC = '_'.join(NCBI_AC)
    except:
        NCBI_AC = a
        #print "Cleaning yielded an error. Got %s, AC is %s" % (a, NCBI_AC)
    return NCBI_AC

def trick(cell):
    if "GCA" in cell:
        cell = cell.replace("GCA", "GCF")
    elif "GCF" in cell:
        cell = cell.replace("GCF", "GCA")
    else:
        cell = cell
    return cell

def getRelevantRank(taxonomy):
    relevant_rank = ""
    #taxonomy = taxonomy.split(";")
    try:
        if "Proteobacteria" in taxonomy:
            taxonomy = taxonomy.split(";")
            for i in taxonomy:
                if i.startswith("c__"):
                    relevant_rank = i[3:]
        else:
            taxonomy = taxonomy.split(";")
            for i in taxonomy:
                if i.startswith("p__"):
                    relevant_rank = i[3:]
    except:
        print "Taxonomy %s yielded an error" % taxonomy
    if relevant_rank == "":
        print "Taxonomy %s returned an empty relevant rank" % taxonomy
    return relevant_rank

def taxonomic_assignment(kmer_value):
    kmer_value = int(kmer_value)
    path = r'intermediate/jellyfish/viral_mae/k%d' % kmer_value
    allFiles = glob.glob(path + "/*viral_batch*.out")
    frame = pd.DataFrame()
    list_ = []
    for file_ in allFiles:
        df = pd.read_table(file_, header=None, sep = "\t")
        list_.append(df)
    
    frame = pd.concat(list_)
    frame.columns = ['viral_contig.kmer', 'MAE', 'host_prediction']
    
    # read NCBI taxonomic assignments from gtdb-tk data
    NCBI_taxonomy = pd.read_table("data/bac_taxonomy_r86.tsv")
    NCBI_taxonomy.columns = ["AC", "taxonomy"]
    NCBI_taxonomy['AC'] = NCBI_taxonomy['AC'].apply(clean)
    
    # since some AC are shifted from GCA to GCF,
    # create a fake version of the NCBI_taxonomy table
    
    NCBI_taxonomy_fake = NCBI_taxonomy.copy()
    NCBI_taxonomy_fake['AC'] = NCBI_taxonomy_fake['AC'].apply(lambda x: trick(x))
    
    # Some genomes are not present in the gtdb-tk data, we had to do that
    NCBI_unknown_taxonomy_files = glob.glob("intermediate/unknown_genomes/unknown_genomes.unknown_genomes.ls.batch*_gtdbtk_out/*classification_pplacer.tsv")
    # Latest batch of taxonomic assignments
    NCBI_unknown_taxonomy_files = NCBI_unknown_taxonomy_files + glob.glob("intermediate/unknown_genomes_2/unknown_genomes.unknown_genomes.ls.batch*_gtdbtk_out/*classification_pplacer.tsv")
    NCBI_unknown_dataframe = pd.DataFrame()
    ulist_ = []
    for file_ in NCBI_unknown_taxonomy_files:
        udf = pd.read_table(file_, header=None, sep = "\t")
        ulist_.append(udf)
    
    uframe = pd.concat(ulist_)
    uframe.columns = ["AC", "taxonomy"]
    
    # read MAG taxonomic assignments
    MAG_taxonomy_files = glob.glob("../gtdb-tk_test/gtdbtk_output_batch.*/*classification_pplacer.tsv")
    #print MAG_taxonomy_files
    mframe = pd.DataFrame()
    mlist_ = []
    for file_ in MAG_taxonomy_files:
        mdf = pd.read_table(file_, header=None, sep = "\t")
        mlist_.append(mdf)
    
    mframe = pd.concat(mlist_)
    mframe.columns = ["AC", "taxonomy"]
    
    # rbind taxonomy tables
    all_taxonomy = NCBI_taxonomy.append(mframe).append(NCBI_taxonomy_fake).append(uframe)
    
    # join!
    frame_taxonomy = pd.merge(frame, all_taxonomy, left_on = "host_prediction", right_on = "AC", how = "left")
    
    # get phylum or class, in case of Proteobacteria
    frame_taxonomy['taxonomy_short'] = frame_taxonomy['taxonomy'].apply(lambda x: getRelevantRank(x))
    os.system('mkdir -p results/viral_contigs_host_prediction/MAE')
    frame_taxonomy.to_csv(path_or_buf = "results/viral_contigs_host_prediction/MAE/viral_contigs_host_assignment_k%d.181114.tsv" % kmer_value, sep = "\t", index = False)
    #return frame_taxonomy

taxonomic_assignment(4)

```

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