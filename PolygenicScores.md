# Building polygenic scores across populations: methods and pitfalls

Scripts and text by Fernando Racimo and Alba Refoyo-Martinez, based on this paper: https://elifesciences.org/articles/39725

In this tutorial, we'll build polygenic scores using trait-associated SNPs from two different GWAS of the same trait: adult human height. We will compute these scores for a panel containing allele frequency data from 26 populations from around the world. Presumably, given that we're looking at the same populations and the same trait, the polygenic scores computed using the two different GWAS should be very similar to each other. We will check if that is so.

The first study is a GWAS meta-analysis conducted by the Genetic Investigation of Anthropometric Traits (GIANT) consortium, in which multiple small GWAS performed across individuals of European ancestry were analyzed together to obtain higher power to test for SNP-trait associations:  
https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files

The second study is a single GWAS performed by the Neale lab on a very large cohort of self-indentified British individuals from the UK Biobank (UKBB), which have largely homogeneous ancestry:
https://www.nealelab.is/uk-biobank/ukbround2announcement

We will obtain P-values denoting the evidence for a SNP-trait association from these studies, along with effect sizes denoting the magnitude and direction of an allele's association to the trait. To compute polygenic scores for a particular population, we also need the allele frequency of each trait-associated allele in that population. To obtain allele frequency data, we will use The 1000 Genomes Project panel, which includes hundreds of individuals sampled from 26 different populations across the world: http://www.internationalgenome.org/category/population/ 

The relevant GWAS and 1000 Genomes files have already been downloaded to our servers, so there is no need to download them yourself. 

First, we'll begin by defining a folder where we've placed scripts to compute these scores:
```
PIPELINEFOL="/home/fernando/PolyScores"
```

We'll also define another folder to place our output files. Where you put them is up to you, just make sure it's within your working directory:
```
OUTPUTFOL=[ here write down a folder to place your output files]
```

We also need a file defining the names of the 26 populations in The 1000 Genomes Project for which we will build polygenic scores. This file is located here:
```
POPS=$PIPELINEFOL/"pops_to_search.txt"
```

# 1 - Using genomic partitions

When building polygenic scores, we want to make sure that the SNPs we use are not in high LD with each other. Some polygenic score methods use information about LD patterns from all SNPs across the genome, to correct for correlations in allele frequencies that may be due to LD. In our case, we will take a conservative approach: we will partition the genome into very large, approximately independent blocks, and use a single SNP from each block to compute our scores. This will ensure that each SNP we use is not in high LD with any other SNP we use. We have already obtained a file containing this block partitions, and it is located here:

```
LDBFILE=$PIPELINEFOL/"fourier_ls-all_nochr.bed"
```

The LD blocks were obtained using a method described here: https://academic.oup.com/bioinformatics/article/32/2/283/1743626


# 2 - Extracting the candidate SNPS

We'll now proceed to compute scores using P-values and effect size estimates from the GIANT study. Let's create a folder to place these results:
```
GWAS="GIANT"
mkdir $OUTPUTFOL/$GWAS
```

The GWAS SNP effect sizes and P-values have already been merged with the 1000 Genomes allele frequencies into a single file. This file cna be found here:
```
RAWGWASFREQ=$PIPELINEFOL/$GWAS/"gwasfreqs_height.tsv"
```
Take a look at thils file (use 'less' in the unix command line). The first 10 columns of this file contain important information about each SNP:

```
zless $PIPELINEFOL/$GWAS/"gwasfreqs_height.tsv.gz"
```

CHROM = chromosome

POS = position (in bp) along the chromosome

SNPID = a unique ID for each SNP

REF = reference allele

ALT = alternative (non-reference) allele

ANC = inferred ancestral (chimp-like) allele

DER = inferred derived (mutant) allele

DEREFFECT = estimated size of the effect of the SNP on the trait (in our case, height), polarized with respect to the derived allele. A positive value denotes that the derived allele is associated with an increase in the trait. A negative value denotes that it is associated with a decrease in the trait.

SE = standard error of the effect size estimate

PVAL = P-value denoting the evidence in favor of an association between the SNP and the trait

The rest of the columns contain allele frequency information for each of the population panels in the 1000 Genomes Project. They are represented as two numbers separated by a comma. The first number is the number of ancestral alleles in that population, while the second number is the number of derived alleles.

From this file, we will extract the SNP with the lowest p-value from each approximately-independent block of the genome. In the script below, the option -p serves to define the maximum p-value cutoff that is allowed for each SNP. If a block does not have a SNP with a P-value lower than the one specified with this option, then it will be ignored. SNPs with smaller p-values are better association candidates, but we also need a high number of SNPs for the polygenic score to be reasonably predictive. In this case, we will use the standard genome-wide P-value significance cutoff: 5e-8. The option -i denotes the input SNP file, the option -b denotes the block file and the option -o denotes the ouput file, which will contain one SNP per block:
```
CANGWASFREQ=$OUTPUTFOL/$GWAS/"gwasfreqs_candidates_height.tsv"
python $PIPELINEFOL/partitionUKB_byP.py -i $RAWGWASFREQ.gz -b $LDBFILE -o $CANGWASFREQ -p5e-08
```

Take a look at your output folder. It should now contain a new file with the candidate SNPs. How many SNPs did you extract?

# 3 - Computing the scores

Ok, now we are ready to calculate the polygenic scores:
```
GENSCORES=$OUTPUTFOL/$GWAS/"Genscores_height.txt"
Rscript $PIPELINEFOL/PolygenicScores.R -w $CANGWASFREQ -p $POPS -s $GENSCORES
```
Here, the option -w denotes the candidate SNP file, the option -p denotes the file containing the list of populations for which we will compute polygenic scores and the option -s denotes the output file we are producing, which contains the polygenic scores for each population. 

You should now have a new file in your output folder, containing the names of each population in the first column and the polygenic scores for height in the second column. What do you observe? Which populations have high polygenic scores? Which ones have low scores? It will help to visualize the results if you plot the scores in R. You can check out the 1000 Genomes Project link at the top of this page, to see a description of what each population label means.

Now, try repeating this exercise but using SNPs determined to be significant in the UK Biobank GWAS instead, and using effect sizes derived from this GWAS as well. To do so, just repeat Step 2 and Step 3, but this time, change the name of the GWAS that you will use:

```
GWAS="UKBB"
```

# 4 - Conclusions

Try plotting the two sets of polygenic scores (obtained from the GIANT and UK Biobank GWAS) on the same plot, using R. What do you observe? Do you observe differences in scores between populations? Are these differences consistent across the two GWAS? Do the differences have the same magnitude regardless of the GWAS used? What could be the cause for inconsistencies between the two approaches?

For further information, you can read the papers that first described these inconsistencies here:

https://cdn.elifesciences.org/articles/39725/elife-39725-v1.pdf

https://cdn.elifesciences.org/articles/39702/elife-39702-v1.pdf

A summary of their findings can be found here: 

https://elifesciences.org/articles/45380

