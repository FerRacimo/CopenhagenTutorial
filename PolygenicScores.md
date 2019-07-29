Scripts and text by Fernando Racimo and Alba Refoyo-Martinez

PIPELINEFOL="/home/fernando/PolyScores"

GWAS="GIANT"

# Data and output folders - modify these to set them to your working environment
# Input file
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
