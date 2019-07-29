# Building polygenic scores across populations: methods and pitfalls

Scripts and text by Fernando Racimo and Alba Refoyo-Martinez

In this tutorial, we'll build polygenic scores using trait-associated SNPs from two different GWAS studies of the same trait: adult human height. We will compute these scores for the population panel, containing allele frequency data from 26 populations from around the world. Presumably, given that we're looking at the same populations and the same trait, the polygenic scores computed using the two different GWAS should be very similar to each. We will check if that is the case.

The first GWAS study is a meta-analysis conducted by the Genetic Investigation of Anthropometric Traits (GIANT) consortium, in which multiple small GWAS performed across individuals of European ancestry were analyzed together to obtain higher power to test for SNP-trait associations:  
https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files

The second GWAS study is a single GWAS performed by the Neale lab on a very large cohort of self-indentified British individuals from the UK Biobank, which have largely homogeneous ancestry:
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

We also need a file defining the names of the populations in The 1000 Genomes Project for which we will build polygenic scores. This file is located here:
```
POPS=$PIPELINEFOL/"pops_to_search.txt"
```

When building polygenic scores, we also want to make sure that the SNPs we use are not in high LD with each other. Some polygenic score methods use information about LD patterns across all SNPs across the genome, to correct for correlations in allele frequencies that may be due to LD. In our case, we will take a conservative approach: we will partition the genome into very large, approximately independent blocks, and use a single SNP from each block to compute our scores. This will ensure that each SNP we use is not in high LD with any other SNP we use. We have already obtained a file containing this block partitions, and it is located here:

```
LDBFILE=$PIPELINEFOL/"fourier_ls-all_nochr.bed"
```


Output file names
```
GWAS="GIANT"
RAWGWASFREQ=$PIPELINEFOL/$GWAS/"gwasfreqs_height.tsv"
CANGWASFREQ=$PIPELINEFOL/$GWAS/"gwasfreqs_candidates_height.tsv"
GENSCORES=$PIPELINEFOL/$GWAS/"Genscores_height.txt"
```

We will first extract the SNP with the lowest p-value from each LD block. The option -p serves to define the maximum p-value cutoff that is allowed for each SNP. If a block does not have a SNP with a P-value lower than the one specified with this option, then it will be ignored. SNPs with smaller p-values are better association candidates, but we also need a good number of SNPs to compute the polygenic score. In this case, we will use the standard genome-wide P-value significance cutoff: 5e-8.
```
python $PIPELINEFOL/partitionUKB_byP.py -i $RAWGWASFREQ.gz -b $LDBFILE -o $CANGWASFREQ -p5e-08
```

Ok, now we are to calculate the polygenic scores:
```
Rscript $PIPELINEFOL/PolygenicScores.R -w $CANGWASFREQ -p $POPS -s $GENSCORES
```

Try repeating this exercise but using SNPs determined to be signifcant in the UK Biobank GWAS instead, and using effect sizes derived from this GWAS as well.

Set:
```
GWAS="UKB"
```
