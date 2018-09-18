#!/bin/sh
# Prints evidence for each base from alignment
# Functions

function msg
{
	echo -e "$*"
}

function err
{
	echo "ERROR: $*" 1>&2
	exit 1
}

function usage
{
	msg "Name:\n  extract_base_evidence"
	msg "Description:\n  Prints evidence for each base from alignment"
	msg "Author:\n  Jason Kwong <jason.kwong@unimelb.edu.au>"
	msg "Usage:\n  extract_base_evidence [options] snps.bam > evidence.txt"
	msg "  Output displayed in tab-delimited format:"
	msg "  CHR	POS	A	C	G	T"
	msg "Parameters:"
	msg "  snps.bam   Sorted alignment"
	msg "Options:"
	msg "  -h         Show this help"
	exit 0
}

#..............................................................................

# Options
while getopts "h" opt; do
	case $opt in
		h)
			usage ;;
	esac
done

# skip over options to pass arguments
shift $((OPTIND - 1))

# read our mandatory positional parameters
[[ $# != 1 ]] && usage

# check whether input file exists
[[ ! -f $1 ]] && err "Input alignment file not found."

samtools mpileup --no-BAQ --min-MQ 60 --min-BQ 13 $1 | cut -f 1,2,5 | perl -ne '@x=split m/\t/; $A = $x[2] =~ tr/Aa/Aa/; $C = $x[2] =~ tr/Cc/Cc/; $G =
$x[2] =~ tr/Gg/Gg/; $T = $x[2] =~ tr/Tt/Tt/;  print join("\t", $x[0],$x[1],$A,$C,$G,$T),"\n";'
