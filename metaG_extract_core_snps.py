#!/usr/bin/env python3
# Extracts most frequent base call from metagenome pileup evidence to add to core.aln

# Usage
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Extracts most likely base call from pileup evidence to add to core.aln\n'
		'Use with "extract_base_evidence.sh"',
	usage='\n  %(prog)s [--out invariant.fa] EVIDENCE')
parser.add_argument('evidence', metavar='EVIDENCE', nargs=1, help='output file from "extract_base_evidence"')
parser.add_argument('--id', metavar='ID', required=True, help='sample/strain ID')
parser.add_argument('--sites', metavar='FILE', required=True, help='file with list of core SNP sites; can use core.tab or core.vcf')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
args = parser.parse_args()

sites_ext = os.path.splitext(args.sites)[-1].lower()

chr = []
sites = []
seq = []
bases = {0:'A', 1:'C', 2:'G', 3:'T'}

if sites_ext == ".tab" or ".vcf":
	with open(args.sites, 'r') as f:
		while True:
			line = f.readline()
			if not line.startswith("#") and not line.startswith("CHR"):
				chrpos = line.split('\t')
				chr.append(chrpos[0])
				sites.append(chrpos[1])
				break
		for line in f.readlines():
			chrpos = line.split('\t')
			chr.append(chrpos[0])
			sites.append(chrpos[1])

core_sites = pd.DataFrame({'locus':chr, 'site':sites})
chr_locus = core_sites['locus'][0]
chr_sites = core_sites[core_sites['locus'] == chr_locus]
core_sites_site = chr_sites['site'].tolist()
#print(chr_sites)

df = pd.read_csv(args.evidence[0], delimiter='\t', header=None, names=['locus', 'site', 'A', 'C', 'G', 'T'])
df2 = df[df['site'].isin(chr_sites['site'])]
df2 = df2[df2['locus'] == chr_locus]
temp_sites = df2['site'].tolist()

# Check missing sites with 0 depth
missing = []
for i in core_sites_site:
	if int(i) not in temp_sites:
		missing.append(i)
temp_locus = [chr_locus] * len(missing)
empty_list = [int(0)] * len(missing)
col_names = ['locus', 'site', 'A', 'C', 'G', 'T']
missing_sites = pd.DataFrame({'locus':temp_locus, 'site':missing, 'A':empty_list, 'C':empty_list, 'G':empty_list, 'T':empty_list})
missing_sites = missing_sites[col_names]

# Append missing_sites df to core_sites df
df3 = df2.append(missing_sites, ignore_index=True)
df3[['site', 'A', 'C', 'G', 'T']] = df3[['site', 'A', 'C', 'G', 'T']].apply(pd.to_numeric)
df3 = df3.sort_values(by=['site'])
df3 = df3.reset_index(drop=True)
#print(df3)

# Identify most abundant base call
for index, row in df3.iterrows():
	acgt = [row['A'], row['C'], row['G'], row['T']]
	if sum(acgt) == 0:
		base = '-'
	else:
		ix = acgt.index(max(acgt))
		base = bases[ix]
	seq.append(base)

# Print sequence to stdout
print(">" + args.id)
print(''.join(seq))
