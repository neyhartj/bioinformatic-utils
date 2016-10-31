# bioinformatic-utils
Scripts to perform various bioinformatics operations. These are generally
written in Python, but I use whatever language is most apropriate.

### Conversion Tools

#### fastPHASE

`tassel2fastPHASE.py` converts a TASSEL hapmap file to a fastPHASE input file.
`fastPHASE2tassel.py` converts a fastPHASE output file to a TASSEL hapmap file. It
requires the original TASSEL hapmap file used to convert to fastPHASE input.


`fastPHASE2vcf.py` converts the output of fastPHASE to a VCF. It requires
the output from fastPHASE (i.e. "switch_guess") and an original VCF file that
was used to convert to fastPHASE input (PLINK has a setting to convert a VCF
to fastPHASE input). If the low genotype likelihood setting was used in 
fastPHASE, those genotypes will be set to missing.

#### Other

`vcf2hapmap.py` converts a VCF to a hapmap file that is either suitable for 
TASSEL or suitable for the `rrBLUP` R package.


### Filtering and Prep

`LB_impute_vcf_prep.py` filters and prepared VCF files for used in the imputation
program [LB-Impute](https://github.com/dellaporta-laboratory/LB-Impute). It requires
a VCF and a pedigree file. The program first filters parents based on markers that
are homozygous within parents and polymorphic between parents. The program then
evaluates progeny for degree of missingness.
