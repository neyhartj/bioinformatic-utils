#!/usr/bin/env python

## This python script will convert a TASSEL hapmap to a file input 
# used by fastPHASE

# Import libraries
import argparse

#####
# Define the arguments
#####

# Description
DESC = """A Python program to convert a TASSEL-encoded hapmap file to the input
file used by fastPHASE.\n"""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
                '--filein',
                metavar = 'FILEIN',
                help = 'The TASSEL-encoded hapmap file to be converted.',
                required = True)
# Output file name
parser.add_argument('-o',
                '--fileout',
                metavar = 'FILEOUT',
                help = 'The name of the output file. Default extension is ".imp."',
                required = True)


# Parse the arguments
args = parser.parse_args()

# Define the output filename
outfile = str(args.fileout) + '.imp'

# # Open a handle for writing
handle = open(outfile, 'w')

# Read in the TASSEL file
with open(args.filein, 'r') as hmp:

	# Create an empty list to store positions
	snp_pos = ['P']
	# Empty counter for number of snps
	n_snps = 0
	# Empty variable for number of samples 
	n_samples = ""

	# Empty list to store matrix components
	geno_matrix = []

	# Iterate over lines
	for line in hmp:

		# The first line has the line names and we need to extract those
		if line.startswith('rs#'):

			# Chomp and split by tab
			tmp = line.strip().split('\t')

			# Find the position of the first sample name
			sample_start = tmp.index('QCcode') + 1
			# Extract the line names
			names = tmp[sample_start:]

			# Find the length of the names list. This is the number of samples
			n_samples = len(names)

		else:

			# Add to the snp counter
			n_snps += 1

			# Strip and split
			tmp = line.strip().split('\t')

			# Append the snp position to the list
			snp_pos.append(tmp[3])

			# The genotypes occur after the 11th column (10 in base 0)
			genotypes = tmp[11:]
			# Add these to the matrix
			geno_matrix.append(genotypes)

	# Print the number of samples
	handle.write(str(n_samples) + '\n')
	# Print the number of snps
	handle.write(str(n_snps) + '\n')
	# Print the position of snps
	handle.write(" ".join(snp_pos) + '\n')

	# Iterate and count over the names of samples
	for j, name in enumerate(names):

		# First write the sample name to the file
		handle.write("# " + name + '\n')

		# Empty list for storing genotypes from a sample
		sample_genos_pos = []
		sample_genos_neg = []

		# Now iterate across the matrix that we created
		for i in range(n_snps):
	
			geno = geno_matrix[i][j]

			# Split
			geno0 = geno[0]
			geno1 = geno[1]

			# Append or convert
			if geno0 == "N":
				sample_genos_pos.append("?")
			else:
				sample_genos_pos.append(geno0)

			if geno1 == "N":
				sample_genos_neg.append("?")
			else:
				sample_genos_neg.append(geno1)

		handle.write("".join(sample_genos_pos) + '\n')
		handle.write("".join(sample_genos_neg) + '\n')

# # Close the file
handle.close()














			




