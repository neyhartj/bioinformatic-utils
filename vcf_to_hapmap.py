#!/usr/bin/env python

# This python script will convert a VCF file to hapmap files suitable for use in TASSEL
## or rrBLUP. The TASSEL version encodes genotypes by their nucleotides (i.e. AA, AG, or GG), 
## while the rrBLUP version encodes homozygotes for the reference allele as 1, heterozygotes at 0, and
## homozygous alternates as -1

# Import modules
import argparse # To get the arguments
import sys

#####
# Define functions
#####

def rrBLUP_hapmap(VCF):
  print "Writing the non-reformatted hapmap file using rrBLUP encoding."
      # Create filename
  filename = str(args.outfile) + '_hmp.txt'
  # Open handle for writing file
  handle = open(filename, 'w')
      
      # Lists for handling chromosome name
  chrom_l = [] # Empty list of chromosomes
  chrom_index = [] # Index of the chromosome positions

      # Start reading through the vcf file
  for line in VCF:

    if line.startswith('##'):
      continue
      # Look for the #CHROM line as this is the line that contains sample information
    elif line.startswith('#CHROM'):
      # Split the line on tabs
      tmp = line.strip().split('\t')
      # Fine the format field
      format_field = tmp.index('FORMAT')
      # Get the samples out of the list
      # Add 1 to the index, so "FORMAT" isn't included
      samples = tmp[format_field + 1:]

      # First the header line with the information
      handle.write('rs#\tallele\tchrom\tpos\t' + '\t'.join(samples) + '\n')

    # Now we have the sample names; let's move on to the genotype data
    else:
      # Create a new list to store the data
      toprint = []

      tmp = line.strip().split('\t')

      # Assigning variable
      chrom = tmp[0]
      position = tmp[1]
      ref_allele = tmp[3]
      alt_allele = tmp[4]

      # Handling chromosome names
      if chrom in chrom_l:
        pass
      else:
        chrom_l.append(chrom)

      # Assign the index of the chromosome within the unique list of chromosomes
      ## as the name of that chromosome
      chrom_name = str(chrom_l.index(chrom) + 1)
      
      # The genotype data
      genotypes = tmp[9:]

      # Create variable for the output file
      # Create the alleles variable
      alleles = ref_allele + '/' + alt_allele
      # The position variable was already created
      # Create a name for the SNP
      snp_id = 'S' + chrom_name + '_' + position

      # Append to the list each snp_id, alleles, etc
      toprint.append(snp_id)
      toprint.append(alleles)
      toprint.append(chrom_name)
      toprint.append(position)

      for g in genotypes:
        # The genotype string is separated by :
            # The first element of the genotype string is the genotype call
        call = g.split(':')[0]
        # Genotypes are listed as allele1/allele2
        # Assume the genotypes are unphased
        # 0 = ref, 1 = alt 1
            # 0/0 = homo ref, 0/1 = het, 1/1 = homo alt
                                
      # Encode genotypes in rrBLUP format
        individual_call =''

        if call == '0/0': # If the call is 0/0, declare as 1
          individual_call += '1'
        elif call == '0/1': # If the call is 0/1, declare as 0
          individual_call += '0'
        elif call == '1/1': # If the call is 1/1, declare as -1
          individual_call += '-1'
        else:
          individual_call += 'NA' # If it isn't any of the above, it its missng
        
        # Append the individual calls to the genotype matrix row
        toprint.append(individual_call)

      # Print the organized list
      handle.write('\t'.join(toprint) + '\n')

  print "File was written as " + filename
  # Close the handle
  handle.close()
##### End of function #####

def TASSEL_hapmap(VCF):
  print "Writing the non-reformatted hapmap file using TASSEL encoding."
  # Create filename
  filename = str(args.outfile) + '_hmp.txt'
  # Open handle for writing file
  handle = open(filename, 'w')
  
  # Lists for handling chromosome name
  chrom_l = [] # Empty list of chromosomes
  chrom_index = [] # Index of the chromosome positions

  # Start reading through the vcf file
  for line in VCF:

    if line.startswith('##'):
      continue
      # Look for the #CHROM line as this is the line that contains sample information
    elif line.startswith('#CHROM'):
      # Split the line on tabs
      tmp = line.strip().split('\t')
      # Fine the format field
      format_field = tmp.index('FORMAT')
      # Get the samples out of the list
      # Add 1 to the index, so "FORMAT" isn't included
      samples = tmp[format_field + 1:]

      # First the header line with the information
      handle.write('rs#\tallele\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode' + '\t'.join(samples) + '\n')

    # Now we have the sample names; let's move on to the genotype data
    else:
        # Create a new list to store the data
      toprint = []

      tmp = line.strip().split('\t')

      # Assigning variable
      chrom = tmp[0]
      position = tmp[1]
      ref_allele = tmp[3]
      alt_allele = tmp[4]

      # Handling chromosome names
      if chrom in chrom_l:
        pass
      else:
        chrom_l.append(chrom)

      # Assign the index of the chromosome within the unique list of chromosomes
      ## as the name of that chromosome
      chrom_name = str(chrom_l.index(chrom) + 1)
      
      # The genotype data
      genotypes = tmp[9:]

      # Create variable for the output file
      # Create the alleles variable
      alleles = ref_allele + '/' + alt_allele
      # The position variable was already created
      # Create a name for the SNP
      snp_id = 'S' + chrom_name + '_' + position

      # Append to the list each snp_id, alleles, etc
      toprint.append(snp_id)
      toprint.append(alleles)
      toprint.append(chrom_name)
      toprint.append(position)
      toprint.extend(['NA'] * 7)

      for g in genotypes:
        # The genotype string is separated by :
        # The first element of the genotype string is the genotype call
        call = g.split(':')[0]
        # Genotypes are listed as allele1/allele2
        # Assume the genotypes are unphased
        # 0 = ref, 1 = alt 1
        # 0/0 = homo ref, 0/1 = het, 1/1 = homo alt
                                
        # Encode genotypes in rrBLUP format
        individual_call =''

        if call == '0/0': # If the call is 0/0, declare as 1
          individual_call += ref_allele + ref_allele
        elif call == '0/1': # If the call is 0/1, declare as 0
          individual_call += ref_allele + alt_allele
        elif call == '1/1': # If the call is 1/1, declare as -1
          individual_call += alt_allele + alt_allele
        else:
          individual_call += 'NN' # If it isn't any of the above, it its missng
          
        # Append the individual calls to the genotype matrix row
        toprint.append(individual_call)

      # Print the organized list
      handle.write('\t'.join(toprint) + '\n')
    
  print "File was written as " + filename
  # Close the handle
  handle.close()
##### End of function #####


#####
# Define the arguments
#####

# Description
DESC = """A Python program to convert a VCF file to a hapmap file in rrBLUP format {-1, 0, 1}. 
This tool is part ofthe GBarleyS pipeline.\n\n"""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
  '--vcf_in',
  metavar = 'VCFIN',
  help = 'Input VCF file',
  required = True)
# Output file name
parser.add_argument('-o',
  '--outfile',
  metavar = 'OUTFILE',
  help = 'Output file basename (i.e. no extension)',
  required = True)
group = parser.add_mutually_exclusive_group(required = True)
group.add_argument('-r',
  '--rrBLUP',
  action = 'store_true',
  help = 'Boolean flag for whether a hapmap file should be exported in rrBLUP format',
  required = False)
group.add_argument('-t',
  '--TASSEL',
  action = 'store_true',
  help = 'Boolean flag for whether a hapmap file should be exported in TASSEL format',
  required = False)

# Parse the arguments
args = parser.parse_args()


#####
# Execute program functions
#####

# Convert to rrBLUP
if args.rrBLUP:
    with open(args.vcf_in, 'r') as VCF:
        rrBLUP_hapmap(VCF)

elif args.TASSEL:
    with open(args.vcf_in, 'r') as VCF:
        TASSEL_hapmap(VCF)
else:
    print "\nNo hapmap arguments given.\n"
    parser.print_help()
    exit(1)
