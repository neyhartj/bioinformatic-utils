#!/usr/bin/env python

## This python script will convert a TASSEL hapmap to a file input 
# used by fastPHASE

# Import libraries
import argparse

#####
# Define the arguments
#####

# Description
DESC = """A Python program to convert the output from fastPHASE to a TASSEL-encoded
hapmap file.\n"""

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
parser.add_argument('-t',
				'--tinput',
				metavar = 'TIN',
				help = 'The original TASSEL file converted to the fastPHASE input file.',
				required = True)
# Output file name
parser.add_argument('-o',
                '--fileout',
                metavar = 'FILEOUT',
                help = 'The name of the output file. Default extension is ".hmp.txt',
                required = True)


# Parse the arguments
args = parser.parse_args()

# Define the output filename
outfile = str(args.fileout) + '.hmp.txt'

# # Open a handle for writing
handle = open(outfile, 'w')

# Empty list to store the initial TASSEL data
snp_info = []

# Read in the original fp input file to obtain snp positions
with open(args.tinput, 'r') as tin:

	# Iterate over lines
	for line in tin:

		# Skip the first line
		if line.startswith('rs#'):
			continue

		else:

			# We need the SNP name, allele, chromosome, and position
			tmp = line.strip().split('\t')

			# Append this information
			snp_info.append(tmp[0:4])

# Find the number of snps
n_snps = len(snp_info)


# List to store column names / sample names
columns = ['rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#', 'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode']

# Empty list to store matrix components
geno_matrix = []

s

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

		# If the line switch is one, begin parsing
		else:

			# If the line begins with #, it is a sample name
			if line.startswith("#"):

				# Strip and split on spaces
				tmp = line.strip().split(' ')

				# The sample name is the second entry
				name = tmp[1]
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


# Print the column names line
handle.write('\t'.join(columns) + '\n')


## Now lets play around with this new geno_matrix
# Iterate over the number of snps
for i in range(len(snp_info)):

	# Empty toprint list
	toprint = []

	# Add the snp_info
	toprint.extend(snp_info[i])

	# Buffer with NAs
	toprint.extend(['NA'] * 7)

	# Collect the genotype information for the i-th snp across all samples
	# Use slice notation to go every other entry
	for j in range(len(geno_matrix))[::2]:

		# Collect allele information
		allele1 = geno_matrix[j][i]
		allele2 = geno_matrix[j+1][i]

		# Convert to N if necessary
		if allele1 == "?":
			allele1 = "N"
		if allele2 == "?":
			allele1 = "N"

		# Concatenate allele information
		geno = str(allele1) + str(allele2)

		# Extend the toprint list
		toprint.append(geno)

	# Write the list to the handle
	# print '\t'.join(toprint)
	handle.write('\t'.join(toprint) + '\n')

# Close the file
handle.close()