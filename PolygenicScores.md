# Introduction

Scripts and text by Fernando Racimo and Alba Refoyo-Martinez

In this tutorial, we'll build polygenic scores using trait-associated SNPs from two different GWAS studies of the same trait: adult human height. The first study is a meta-analysis conducted by the Genetic Investigation of Anthropometric Traits (GIANT) consortium, in which multiple small GWAS were analyzed together to obtain higher power to test for SNP-trait associations. The second study is a single GWAS performed on a very large cohort of self-indentified British individuals from the UK Biobank, which have largely homogeneous ancestry. 

First, we'll begin by defining a folder where we've placed scripts to compute these scores
PIPELINEFOL="/home/fernando/PolyScores"

We'll also define another folder to place our output files. Where you put them is up to you, just make sure it's within your working directory.
OUTPUTFOL=[ here write down an folder to place your output files]

GWAS="GIANT"

LDBFILE=$PIPELINEFOL/"fourier_ls-all_nochr.bed"
POPS=$PIPELINEFOL/"pops_to_search.txt"

# Output file names
RAWGWASFREQ=$PIPELINEFOL/$GWAS/"gwasfreqs_height.tsv"
SELGWASFREQ=$PIPELINEFOL/$GWAS/"gwasfreqs_candidates_height.tsv"
GENSCORES=$PIPELINEFOL/$GWAS/"Genscores_height.txt"

# Populations to study (comma-separated) - you can include as many as you want
POPS=$PIPELINEFOL/"pops_to_search.txt"

# Extract best SNP per LD block - maximum p-value cutoff is optional (SNPs with smaller p-values are better association candidates, but you also want a good number of SNPs)
# max pvalue allowed genome wide significance pval
python $PIPELINEFOL/partitionUKB_byP.py -i $RAWGWASFREQ.gz -b $LDBFILE -o $SELGWASFREQ -p5e-08

# 4 - CALCULATION OF STATISTICS

# Calculate polygenic scores
echo 'Polygenic scores'
Rscript $PIPELINEFOL/PolygenicScores.R -w $SELGWASFREQ -p $POPS -s $GENSCORES
