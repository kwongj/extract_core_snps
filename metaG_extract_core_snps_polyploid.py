#!/usr/bin/env python3
# Extracts most frequent base call from metagenome pileup evidence to add to core.aln

# Usage
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import re
import pandas as pd
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from io import StringIO

# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1)

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Extracts most likely base call from pileup evidence to add to core.aln\n'
		'Use with "extract_base_evidence.sh"',
	usage='\n  %(prog)s [--out invariant.fa] EVIDENCE > id.core.snps')
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

msg('Reading base evidence from {} ...'.format(args.evidence[0]))
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

dfcore = pd.DataFrame({'locus':chr, 'site':sites})
all_loci = dfcore['locus'].tolist()
n_loci = len(all_loci)
set_loci = list(OrderedDict.fromkeys(all_loci))

core_loci = []
for index, row in dfcore.iterrows():
	chr_site = str(row[0]) + ':' + str(row[1])
	core_loci.append(chr_site)

df = pd.read_csv(args.evidence[0], delimiter='\t', header=None, names=['locus', 'site', 'A', 'C', 'G', 'T'])
df['chr_site'] = df.locus.astype(str).str.cat(df.site.astype(str), sep=':')
df2 = df[df['chr_site'].isin(core_loci)]
cov_loci = df2['chr_site'].tolist()

missing = []
missing_loci = []
missing_sites = []
for chr_site in core_loci:
	msg('Reading core SNP site: ' + str(core_loci.index(chr_site)+1), end='\r')
	if chr_site not in cov_loci:
		locus = chr_site.split(':')[0]
		site = chr_site.split(':')[1]
		missing.append(chr_site)
		missing_loci.append(locus)
		missing_sites.append(site)
#print(len(missing)+len(cov_loci))

empty_list = [int(0)] * len(missing)
dfmissing = pd.DataFrame({'locus':missing_loci, 'site':missing_sites, 'A':empty_list, 'C':empty_list, 'G':empty_list, 'T':empty_list, 'chr_site':missing}, columns=['locus', 'site', 'A', 'C', 'G', 'T', 'chr_site'])
#print(dfmissing)

# Append missing_sites df to core_sites df
df3 = df2.append(dfmissing, ignore_index=True)
df3[['site', 'A', 'C', 'G', 'T']] = df3[['site', 'A', 'C', 'G', 'T']].apply(pd.to_numeric)
df3 = df3.sort_values(by=['site'])
df3 = df3.reset_index(drop=True)
msg('\nTotal core sites: ' + str(len(df3['chr_site'].tolist())))

# Identify most abundant base call
for index, row in df3.iterrows():
	acgt = [row['A'], row['C'], row['G'], row['T']]
	if sum(acgt) == 0:
		base = '-'
	else:
		ix = acgt.index(max(acgt))
		base = bases[ix]
	seq.append(base)
newseq = ''.join(seq)
newseq = Seq(re.sub('[^ACGTN-]',"",newseq.upper()))

# Print sequence to stdout
seqWRITE = SeqRecord(newseq, id=args.id, description='')
seqFILE = StringIO()
SeqIO.write(seqWRITE, seqFILE, 'fasta')
output = seqFILE.getvalue().rstrip()
print(output)

