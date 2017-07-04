# TE_HUNTER
A pipeline to find TE insertions from NGS data
Version alpha 1.0-2017-07

-----------------
## Author : 
DECHAUD Corentin

@E-mail : corentin.dechaud@gmail.com

-----------------
## Dependencies :
  - Python2.7 :
You need python2.7 to run those scripts. If you don't have it you can install it by running :
```
sudo apt-get install python
```
  - Argparse module for python3 :
 You can [download](https://pypi.python.org/pypi/argparse#downloads) argparse here.
 
  - re module for python2.7.
  
  - Python libraries : sys, os, subprocess, string and threading used in the script are delivered with python. 
  
  - Blast : blastall need to be accessed in the PATH.
  
  - Bowtie2 : Needs to be installed in `/usr/remote/bin/`. In the future some option may be created to seek it in another directory.
  
-----------------
## Install :
First clone this repository on your computer.

```
mkdir /my/path/to/repository
cd /my/path/to/repository
git clone https://github.com/Co-Dec/TE_HUNTER
cd TE_HUNTER
chmod +x pipeline_tehunter.py
cd
```

------------------------
## General description :

Based on different files presented after, this pipeline finds TE insertions in a genome. It can find TE annotated in a reference genome or de novo insertions that are not annotated.

Input data :
  - Single end NGS data.
  - Reference genome.
  - Reference with TE masked genome.
  - TE annotation in the reference.

#### First step :

Retrieve the split-reads, mapping on TEs and genome.

#### Second step :

Split the reads, by keeping only their genomic parts.

#### Third step :

Assemble those splitted reads.

#### Fourth step :

Map those contigs on the reference genome.

#### Fifth step *TE_Hunter.png* :

Check if those contigs are close to annotated TEs. Or if 2 contigs are close without a TE between them.

-----------------
## How to run the script :
Go in the directory where you want to work, could be different directory than the script.

```
cd my/working/dir/
```

In this directory you need madatory files :

  - Single end NGS data : (-reads)
  
The file needs to be cleaned before the analysis. FASTQ format is madatory. I recommend a coverage between 30X and 80X.

  - Reference genome : (-ref)
  
The reference genome for the species studied. FASTA format is mandatory. The header line of each chromosome needs to be splitted by spaces, and in the first space the chromosome name is required. For example : 

**>2L type=golden_path_region; loc=2L:1..23513712; ID=2L; dbxref=GB:AE014134,GB:AE014134,REFSEQ:NT_033779; MD5=b6a98b7c676bdaa11ec9521ed15aff2b; length=23513712; release=r6.16; species=Dmel;**

  - Reference masked genome : (-mask)
  
The reference genome for the studied species but with TE masked and replaced by "N". We recommend to use either *RepeatMasker* or a script present in this package : `TE_masker.py`. It uses the reference genome and the TE sequences file in flybase format.

If you choose to use `RepeatMasker`, please only mask TE, don't mask all repeated regions.

  - TE annotation : (-tebank)
  
The TE annotation consists of the sequences of the TE annotated in the reference genome.

1 : You have those TE in a FASTA format with headers from flybase in this format :

**
>FBti0019256 type=transposable_element; loc=2L:22300300..22304444; name=invader2{}555 dbxref=FlyBase_Annotation_IDs:TE19256,FlyBase:FBti0019256; MD5=d9259a0e33aad699215e64916bd47a5b; length=4145; release=r6.16; species=Dmel;**

Here you can use this file directly and start the pipeline.

2 : You don't have this kind of headers. So you need to change your headers so they look like this :

**>FBtiNUMBER_HERE type=transposable_element; loc=CHR_HERE:0...0; name=FAMILY_HERE{}; dbxref=_; MD5=NoMD5; length=LENGTH_HERE;**

You need to fill NUMBER / CHR / FAMILY and LENGTH. Either you provide TE position in the "loc" field, either you have to add another file : (-ted)

`TE_Description.csv` : A CSV file, separated by tabulations, with : FBti / Family / Chromosome / Start / End.
If you don't provide it and put the location in the headers of TE sequences, this file will be generated by the pipeline.

  - Working dir : (-workdir)
  
The place where you can find your input files, exemple : `/my/working/dir/`

  - Prefix : (-preftoremove)
  
A prefix for all of the files created by the pipeline so you can remove them easily, exemple : `Run1MySpecies30X`

Once all of these files are created you can run TE_HUNTER :

  - Cores : (-cores)
  
Number of threads you want to use.

```
/my/path/to/repository/TE_HUNTER/pipeline_tehunter.py -ref dmel-all-chromosome-r6.16.fasta -cores 16 -reads LibPE_DmGoth10-1.fastq -tebank dmel-all-transposon-r6.16.fasta -mask dmel-all-masked.fasta -workdir /pandata/dechaud/Dmen_verif/ -preftoremove DmGoth10-1
```

## On a calcul cluster :

Watch out this pipeline take a lot of ram especially when using blast :

```
#!/bin/bash

#PBS -q q1day
#PBS -l nodes=1:ppn=16,mem=12gb
#PBS -e /pandata/dechaud/Dmen_verif/DmGoth.err
#PBS -o /pandata/dechaud/Dmen_verif/DmGoth.out
#PBS -N dmel_goth

hostname

/pandata/dechaud/Ancien/TE_st_test/pipeline_tehunter.py -ref dmel-all-chromosome-r6.16.fasta -cores 16 -reads LibPE_DmGoth10-1.fastq -tebank dmel-all-transposon-r6.16.fasta -mask dmel-all-masked.fasta -workdir /pandata/dechaud/Dmen_verif/ -preftoremove DmGoth10-1
```
 

-----------------
## Contacts :
DECHAUD Corentin

E-mail : corentin.dechaud@gmail.com
