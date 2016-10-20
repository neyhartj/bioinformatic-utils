#!/usr/bin/env python

## This python script will convert a TASSEL hapmap to a file input 
# used by fastPHASE

# Import libraries
import argparse
from datetime import datetime
import sys

#####
# Define the arguments
#####

# Description
DESC = """A Python program to convert the output from fastPHASE to a VCF.\n"""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input file
parser.add_argument('-i',
                '--filein',
                metavar = 'FILEIN',
                help = 'The fastPHASE output file to be converted.',
                required = True)
# The original TASSEL file
parser.add_argument('-v',
				'--vcfin',
				metavar = 'VCFI',
				help = 'The original VCF file of unimputed genotypes (the number of sites in the VCF must match the number of sites in the fastPHASE output).',
				required = True)
# Output file name
parser.add_argument('-o',
                '--fileout',
                metavar = 'FILEOUT',
                help = 'The basename of the output VCF file.',
                required = True)


# Parse the arguments
args = parser.parse_args()

# Define the output filename
outfile = str(args.fileout) + '.vcf'

# # Open a handle for writing
handle = open(outfile, 'w')

# Empty list to store the chromosome, position, ID, and ref/alt states of the sites
snp_info = []

# Read in the original fp input file to obtain snp positio
with open(args.vcfin, 'r') as vcfin:

	# Iterate over lines
	for line in vcfin:

		# Skip the first line
		if line.startswith('#'):
			continue

		else:

			# We need the SNP name, allele, chromosome, and position
			tmp = line.strip().split('\t')

			# Append this information
			snp_info.append(tmp[0:5])

## The snp_info will be used later to print the VCF

## Now we will parse the fastPHASE output file

# List to store column names / sample names
columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

# Empty list to store genotype matrix components (i.e. alleles)
geno_matrix = []

# Read in the fastPHASE file
with open(args.filein, 'r') as fp:

	# Create a switch to activate parsing the genotypic information
	line_switch = False

	# Iterate over lines
	for line in fp:

		if not line_switch:

			# If the line doesn't start with "BEGIN GENOTYPES", skip it
			if line.startswith("BEGIN GENOTYPES"):
				# Flip the switch
				line_switch = True
				continue
			else:
				continue

		# If the line switch is on, begin parsing
		else:

			# If the line begins with #, it is a sample name
			if line.startswith("#"):

				# Strip and split on spaces
				tmp = line.strip().split(' ')

				# The sample name is the third entry
				name = tmp[2]
				# Append it to the columns list
				columns.append(name)

			# The last line contains "END GENOTYPES." Skip it
			elif line.startswith("END GENOTYPES"):
				continue

			# Otherwise the line is genotypic information
			else:
				# Strip and split on spaces
				tmp = line.strip().split(' ')

				# Append to the genotype matrix list
				geno_matrix.append(tmp)


# To format the geno_matrix for a VCF, we need to first iterate
## over entries, then iterate over sites.

# The number of entries is one-half of the length of the geno_matrix (for diploid)
n_entries = len(geno_matrix) / 2
# The number of sites
n_sites = len(snp_info)


# Empty list to store the per-site genotype calls
vcf_matrix = [['.' for i in range(n_entries)] for j in range(n_sites)]

# Iterate over entries
for i in range(n_entries):

	# Create indicies to sample the haplotypes
	ind1 = 2 * i
	ind2 = (2 * i) + 1

	# Pull out the alleles for the ith entry
	hap_i1 = geno_matrix[ind1]
	hap_i2 = geno_matrix[ind2]

	# Iterate over the sites
	for j in range(n_sites):

		# Pull out the site info from the VCF
		site_info = snp_info[j]

		# Pull out the reference / alternate state
		ref = site_info[3]
		alt = site_info[4]

		# Pull out the alleles in the haplotypes of the ith individual
		allele_j1 = hap_i1[j]
		allele_j2 = hap_i2[j]

		# Empty genotype call
		call_j = ''

		# If one of the alleles has a bracket, set to missing
		if '[' in allele_j1 or '[' in allele_j2:

			call_j = './.'

		# Otherwise, match the allele to the reference or alternate state
		else:

			# Match allele 1
			if allele_j1 == ref:
				call_j1 = '0'
			else:
				call_j1 = '1'

			# Match allele 2
			if allele_j2 == ref:
				call_j2 = '0'
			else:
				call_j2 = '1'

			# Combine into the call
			call_j = call_j1 + '|' + call_j2

		# Add the genotype call to the empty vector
		vcf_matrix[j][i] = call_j

## Now we can print out the VCF file
handle = open(outfile, 'w')

# First print some header information
handle.write('##fileformat=VCFv4.2' + '\n')
# Date
date = datetime.now().date().strftime('%Y%m%d')
handle.write('##fileDate=' + date + '\n')
# Others
handle.write('##source=fastPHASEv1.4.8' + '\n')
handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"' + '\n')
handle.write('\t'.join(columns) + '\n')

# Iterate over sites
for j in range(n_sites):

	# Start with the site information
	toprint = snp_info[j]

	# Add other information for QUAL, FILTER, INFO, and FORMAT
	toprint.extend(['.', '.', 'PR', 'GT'])

	# # Add genotype calls
	toprint.extend(vcf_matrix[j])

	# # Write to the handle
	handle.write('\t'.join(toprint) + '\n')
