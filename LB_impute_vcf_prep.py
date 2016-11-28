#!/usr/bin/env python

# This python script will filter a VCF file and prepare it for imputation
## via the software LB-impute. The script does the following actions:
#
# 1. Identify the parent and progeny individuals in a family
# 2. Determine the sites that are heterozygous in either parent - set to missing
# 3. Determine the sites that are not polymorphic between parents - set to missing
# 4. Remove the sites where < y percent of the progeny individuals are non-missing

# Import modules
import argparse # To get the arguments
import sys


####################
# Define Functions #
####################

# Define a function to evaluate the parents
# M: genotype matrix dictionary
# parents: a list of parent names
def parent_eval(M, parents):

	# Notify the user
	print "\nEvaluating sites among the parents..."

	# Create a new dictionary with the updated genotypes (for use in LB-impute)
	correct_mat = {}

	# Create another matrix to keep track of homozygous markers between parents
	homozygous_mat = {}

	# Iterate over sites
	for site in M:

		# At each site, determine a number of things:
		# 1. Are the parents homozygous?
		# 2. Are the parents polymorphic between?
		# 3. Is either parent a het?

		# Create a new dictionary
		parent_genos = dict((parent, M[site][parent]) for parent in parents)

		# Get the genotypes, but maintain the key-value structure
		gt = dict((key, value[0]) for key, value in parent_genos.iteritems())

		# Are both parents the same?
		if len(set(gt.values())) == 1:

			# Are both parents NA?
			if all([call == './.' for call in gt.values()]):

				# Output the same genotypes
				correct_mat[site] = parent_genos

			# Are both parents hets?
			elif all([call == '0/1' for call in gt.values()]):

				# Ouput missing
				new_parent_genos = dict((key, ['./.'] + value[1:]) for key, value in parent_genos.iteritems())
				correct_mat[site] = new_parent_genos

			# The remaining option is that both parents are homozygous for the same allele
			# In this case, output missing, but keep track of these markers
			else:

				# Set both to missing, but keep allele depth
				new_parent_genos = dict((key, ['./.'] + value[1:]) for key, value in parent_genos.iteritems())
				correct_mat[site] = new_parent_genos

				homozygous_mat[site] = parent_genos

		# Otherwise the parents must be different
		else:

			# Is there one missing value?
			if any([call == './.' for call in gt.values()]):

				# Set both to missing, but keep allele depth
				new_parent_genos = dict((key, ['./.'] + value[1:]) for key, value in parent_genos.iteritems())
				correct_mat[site] = new_parent_genos

			# If there isn't a single missing value, both sites must have values (duh)
			# The remaining options here are the following:
			## 0/0 and 0/1 or vice versa (undesirable)
			## 0/0 amd 1/1 or vice versa (desirable)
			## 1/1 and 0/1 or vice versa (undesirable)
			else:

				# Is one genotype a het?
				if any([call == '0/1' for call in gt.values()]):
				
					# Set both to missing, but keep allele depth
					new_parent_genos = dict((key, ['./.'] + value[1:]) for key, value in parent_genos.iteritems())
					correct_mat[site] = new_parent_genos

				# If one is not a het, both genotypes must be homozygous
				else:

					correct_mat[site] = parent_genos

	# Notify
	print "Done\n"

	# Return the corrected matrix
	return correct_mat, homozygous_mat



# Define a function to evaluate the progeny
# M: genotype matrix dictionary
# progenies: a list of progeny names
# max_missing: the maximum missingness proportion on a per-site basis (i.e. across entries)
# This function needs to:
## 1. Find the sites with low missing (i.e. 5 or fewer individuals) and remove
## 2. Find the progeny with > 0.35 heterozygosity
def progeny_eval(M, progenies, max_missing):

	# Notify the user
	print "\nEvaluating sites among the progeny..."

	# Create a new dictionary with the updated genotypes (for use in LB-impute)
	correct_mat = {}

	# Iterate over sites
	for site in M:

		# At each site, determine the following:
		# 1. IS the missingness too high?

		# Create a new dictionary
		progeny_genos = dict((progeny, M[site][progeny]) for progeny in progenies)

		# Get the genotypes
		gt = dict((key, value[0]) for key, value in progeny_genos.iteritems())

		# Count the number of missing
		n_missing = gt.values().count('./.')
		# Total number of gt
		n = len(gt)
		# Proportion of missing
		p_missing = float(n_missing) / n

		# If the p_missing is <= than max_missing, add the genotypes
		if p_missing <= max_missing:

			correct_mat[site] = progeny_genos

	# Notify
	print "Done\n"

	# Return the matrix
	return correct_mat

# End of function

# Define a function to print a new VCF given a dictionary of 
## genotypes and information for each snp
def print_vcf( M, snp_names, snp_info, filename ):

	# Add the extension to the filename
	filename = filename +  '.vcf'

	# Open a file to write
	handle = open(filename, 'w')

	# Write the VCF header
	handle.write('##fileformat=VCFv4.2' + '\n')
	handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n')
	handle.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">' + '\n')
	handle.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">' + '\n')
	
	# Column headers
	headers_toprint = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

	# Extract sample names and add to the headers
	samples = sorted(M.values()[0].keys())
	headers_toprint = headers_toprint + samples

	# Write the VCF column headers
	handle.write('\t'.join(headers_toprint) + '\n')

	# Iterate over sites in the snp_names vector
	for site in snp_names:

		# If the site is in the keys for the matrix, continue
		if site in M.keys():

			toprint = []

			# Find the info in the info_mat and store as a vector to print
			toprint += snp_info[site]

			# Extract the genotype information
			genos_toprint = M[site]

			# Iterate over the sorted keys, extract the values, and add to a vector
			for sample in samples:

				# Extract the genotype
				sample_geno = genos_toprint[sample]

				toprint.append(':'.join(sample_geno))

			# Print
			handle.write('\t'.join(toprint) + '\n')

	# Close the file
	handle.close()

	return

# End of function



######################
# Define the arguments
######################

# Description
DESC = """A Python program to filter and prepare a VCF file for use
in the imputation program 'LB-Impute.' The program requires a VCF file
and a pedigree file that includes parent and progeny IDs. The program 
does two filtering steps. First, it will only keep sites that are homozygous 
within parents, and polymorphic between (all others are set to missing). Next,
the program will filter the progeny for missingness up to a provided threshold.\n\n"""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
  '--vcfin',
  metavar = 'VCFIN',
  help = 'Input VCF file',
  required = True)
# Output file name
parser.add_argument('-o',
  '--vcfout',
  metavar = 'VCFOUT',
  help = 'Output file basename (i.e. no extension)',
  required = True)
parser.add_argument('-p',
  '--pedigree',
  metavar = 'PEDIGREE',
  help = 'Tab-delimited pedigree file. The first column contains individual IDs, the third column contains family IDs, and the third column contains an indicator of "Parent" or "Progeny."',
  required = True)
parser.add_argument('-miss',
  '--missingness',
  metavar = 'MISSINGNESS',
  help = 'The maximum proportion of missingness allowed per-site in the progeny.',
  type = float,
  required = True)

# Parse the arguments
args = parser.parse_args()


#####
# Execute program functions
#####

# Take the VCF from the argument
vcfin = args.vcfin

# Notify the user
print "\nParsing the VCF... "

# Create a matrix of genotype calls
# The matrix will be a dictionary, where the key
# is the site name and the value is the dictionary of entry genotypes
mat = {}
# Create a matrix of SNP information (i.e. pos, id, ref, alt, etc)
info_mat = {}
# Vector of headers
vcf_headers = []
# Vector of column names
headers = []
# Vector of SNP names (keep the order to retrieve later)
snp_names = []

# Empty number of samples scalar
n_entries = 0

# Empty scalar for the index of the first genotypes
entries_index = 0

# Empty vector of entry names
entries = []


# Open and parse the VCF file to extract genotype information to a matrix
with open(vcfin, 'r') as vcf:

	# Iterate over lines
	for line in vcf:

		# Skip the information lines
		if line.startswith("##"):

			vcf_headers.append(line)

		# Extract sample names from the #CHROM row
		elif line.startswith("#"):

			# This is the header vector
			headers = line.strip().split("\t")
			# Extract sample names
			entries_index = headers.index("FORMAT") + 1
			entries = headers[entries_index:]
			# Find the number of entries
			n_entries = len(entries)

		else:

			# Strip and split
			tmp = line.strip().split("\t")

			# Extract SNP information
			info = tmp[:entries_index]

			# Create a name for the SNP using the chromosome and position
			snp_name = info[0] + "_" + str(info[1])
			# Add to the snp_names vector
			snp_names.append(snp_name)

			# Exract or designate other information
			CHROM = info[0]
			POS = info[1]
			ID = '.'
			REF = info[3]
			ALT = info[4]
			QUAL = '.'
			FILTER = '.'
			INFO = '.'
			FORMAT = 'GT:AD:DP'

			# Store #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT
			info_mat[snp_name] = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]

			# Gather the genotypes
			genos = tmp[entries_index:]

			# Create an empty dictionary to store the genotypes for each entry
			site_genos = {}
			
			# Iterate over the entry names and genotype calls
			for entry, geno in zip(entries, genos):

				# Split the genotype information on colon
				call = geno.split(":")

				# Get the genotype call
				gt = call[0]
				# Allele depth
				ad = call[1]
				# Total depth
				dp = call[2]

				# Add the gt and ad to the dictionary
				site_genos[entry] = [gt, ad, dp]

			# Add the site dictionary to the matrix
			mat[snp_name] = site_genos

			

print "Done\n"


# Parse the pedigree file
pedigree = args.pedigree

# Notify
print "\nParsing the pedigree..."

with open(pedigree, 'r') as ped:

	# Create an empty dictionary to store family information
	families = {}

	# Create an empty vector for parennts
	parents = []
	# Create an empty vector for progeny 
	progeny = []

	# Iterate over lines of the pedigree
	for line in ped:

		# Strip and split
		tmp = line.strip().split('\t')

		# First column is individual ID
		individual_id = str(tmp[0])
		# The second column is the family ID
		family_id = str(tmp[1])
		
		# Create a new dictionary for the family, if it doesn't
		## already exit
		if family_id not in families.keys():

			families[family_id] = {'Parents': [], 'Progeny': []}

		# Designate to a vector based on indicator
		if tmp[2] == "Parent":

			# Append the individual ID
			families[family_id]['Parents'].append(individual_id)

			parents.append(individual_id)

		if tmp[2] == "Progeny":

			# Append the individual ID
			families[family_id]['Progeny'].append(individual_id)

			progeny.append(individual_id)

# Notify
print "Done\n"
# Counts
print "There is/are " + str(len(families)) + " family/families, " + str(len(parents)) + " parent(s), and " + str(len(progeny)) + " progeny in the pedigree.\n"


# Separate tasks into families
for family in families:

	# Evaluate the parents
	family_parents = families[family]['Parents']
	eval_parent_genos = parent_eval(mat, family_parents)

	# The first entry in that list is the correct VCF for parents. The second entry
	## are the monomorphic sites (still usefull)
	parent_genos_touse = eval_parent_genos[0]

	#print parent_genos_touse

	parent_genos_mono = eval_parent_genos[1]

	# Subset the progeny
	family_progeny = families[family]['Progeny']
	# Evaluate the progeny
	progeny_genos_touse = progeny_eval(mat, family_progeny, args.missingness)


	# Sites are only removed in the progeny, so we can loop over the sites in the progeny,
	## find the same sites in the parents, and combine
	combined_mat = {}

	for site in progeny_genos_touse:

		# Site dictionary starting with parents
		site_combined = parent_genos_touse[site]
		# Add the progeny
		site_combined.update(progeny_genos_touse[site])

		# Add the site to the larger matrix
		combined_mat[site] = site_combined



	# For some reason the parent monomorphic VCF needs to be written first...
	## Notify
	print "\nWriting the VCF of monomorphic parent genotypes...\n"

	filename = args.vcfout + "_family" + family + "_parent_monomorphic"
	# Print
	print_vcf(parent_genos_mono, snp_names, info_mat, filename)


	# Print the combined VCF
	## Notify
	print "\nWriting the VCF of combined and filtered parent and progeny genotypes...\n"

	## Create a new filename
	filename = args.vcfout + "_family" + family + "_combined"
	# Print
	print_vcf(combined_mat, snp_names, info_mat, filename)


