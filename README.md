# extract_core_snps
Extracts core SNPs from a metagenome for strain-level resolution

## Author

Jason Kwong (@kwongjc)

## Dependencies
* Python 3.x
* Pandas
* Snippy
* bwa mem

## Usage

```
$ extract_base_evidence.sh -h
Name:
  extract_base_evidence
Description:
  Prints evidence for each base from alignment
Author:
  Jason Kwong <jason.kwong@unimelb.edu.au>
Usage:
  extract_base_evidence [options] snps.bam > evidence.txt
  Output displayed in tab-delimited format:
  CHR	POS	A	C	G	T
Parameters:
  snps.bam   Sorted alignment
Options:
  -h         Show this help
```
  
```
$ metaG_extract_core_snps.py -h
usage: 
  metaG_extract_core_snps.py [--out invariant.fa] EVIDENCE

Extracts most likely base call from pileup evidence to add to core.aln
Use with "extract_base_evidence.sh"

positional arguments:
  EVIDENCE      output file from "extract_base_evidence"

optional arguments:
  -h, --help    show this help message and exit
  --id ID       sample/strain ID
  --sites FILE  file with list of core SNP sites; can use core.tab or core.vcf
  --version     show program's version number and exit
```

## Basic instructions

1. Run [`snippy`](https://github.com/tseemann/snippy) and `snippy-core` to identify core genome SNP sites from a group of comparison genomes using specified reference genome

2. Map metagenome to same reference genome using bwa mem
```
bwa index ref.fa
bwa mem -t $CPUS ref.fa ../clipped_R1.fq.gz ../clipped_R2.fq.gz | samtools sort > metaG_sorted.bam
```
3. Extract evidence from pileup for every site where coverage ≥1
```
extract_base_evidence.sh metaG_sorted.bam > metaG.counts
```
4. Extract base calls from core sites in `core.tab` or `core.vcf` from `snippy` output
```
metaG_extract_core_snps.py --id NAME --sites core.tab metaG.counts > metaG.core.snps
```
5. Add metagenome core SNP site calls to core genome SNP alignment
```
cat core.aln metaG.core.snps > new_core.aln
```


## Bugs

Please submit via the [GitHub issues page](https://github.com/kwongj/extract_core_snps/issues).  

## Software Licence

[GPLv3](https://github.com/kwongj/extract_core_snps/blob/master/LICENSE)
