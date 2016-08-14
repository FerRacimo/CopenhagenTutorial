Tutorial of demographic analyses with ngsTools and ANGSD - Copenhagen 2016 course
===============
Adapted from: 
- ngsTools tutorial by Matteo Fumagalli: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
- ANGSD ABBA-BABA tutorial by Thorfinn Korneliussen: http://popgen.dk/angsd/index.php/Abbababa

In this tutorial we will be using ngsTools and ANGSD to compute summary statistics and perform population genetics analyses from low-depth sequencing data.


First, inside your home directory, create a new directory where we will be working and go into it:
```
mkdir TutorialDemo
cd TutorialDemo
```

Data
----------

As an illustration, we will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and a putative ancestral sequence.
BAM files have been downsampled to a mean depth of around 4X.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

The data can be found here:
/ricco/data/fernando/TutorialFiles

Let’s create a shortcut for this directory, as we will refer to it quite a lot in our command lines:
```
DATAFOL=/ricco/data/fernando/TutorialFiles
```

Here, we have 60 BAM files at low depth. The list with BAM files has been written to 'ALL.bamlist’:
```
cat $DATAFOL/ALL.bamlist
```


We also have a reference and an ancestral sequence in FASTA format. Let’s also make shortcuts for them:
```
REF=$DATAFOL/Data/ref.fa.gz
ANC=$DATAFOL/Data/anc.fa.gz
```

As a note for the general use, in case an ancestral sequence is not available, analyses on FST, PCA, nucleotide diversity (but not the number of fixed differences) can be carried out using the reference sequence to polarize your data. Please be aware that, under this scenario, some quantities (e.g. the unfolded joint site frequency spectrum) will be nonsense. Please also note that, since we are randomly subsampling reads here, your results in this tutorial may (slightly) differ from what written here. 

Throughout the tutorial, we will also use a number of R scripts that can be downloaded from the ngsTools github website. Let’s make a shortcut to them too:
```
SCRIPTS=/ricco/data/fernando/Scripts
```

Summary statistics using ANGSD
-----------------------------------

Nucleotide diversity indexes and measures of population differentiation can be estimated taking data uncertainty into account both with ANGSD and ngsTools.
First, we show how to estimate such summary statistics using ANGSD.

## The Site Frequency Spectrum (SFS)

One of the most important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS). 
SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. 
SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described). 
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site. 
Finally, an estimate of the SFS is computed.

These steps can be accomplished in ANGSD using `-doSaf 1/2` options and the program `realSFS`.

```
angsd -doSaf
...
-doSaf		0
	1: perform multisample GL estimation
	2: use an inbreeding version
	3: calculate genotype probabilities
	4: Assume genotype posteriors as input (still beta) 
	-doThetas		0 (calculate thetas)
	-underFlowProtect	0
	-fold			0 (deprecated)
	-anc			(null) (ancestral fasta)
	-noTrans		0 (remove transitions)
	-pest			(null) (prior SFS)
	-isHap			0 (is haploid beta!)
NB:
	  If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site
```

The SFS is typically computed for each population separately.
We need to slightly modify the filtering options as now each population has 20 samples. 
So now we set `-minInd 20 -setMinDepth 20 -setMaxDepth 200`.
Also, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarize our alleles (to ancestral and derived states).

We cycle across all populations:
```
for POP in LWK TSI PEL
do
	echo $POP
	angsd -P 4 -b $DATAFOL/$POP.bamlist -ref $REF -anc $ANC -out $POP \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
		-GL 1 -doSaf 1 &> /dev/null
done
```

Have a look at the output file.
```
realSFS print PEL.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site, as seen during the lecture.
So the first value (after the chromosome and position columns) is the likelihood of having 0 copies of the derived allele, the second indicates the probability of having 1 copy and so on.
Note that these values are in log format and scaled so that the maximum is 0.

The next step would be to use these likelihoods and estimate the overall SFS.
This is achieved by the program `realSFS`.
```
realSFS
	-> ---./realSFS------
	-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:

	-> Estimate the SFS for entire genome??
	-> ./realSFS afile.saf.idx 

	-> 1) Estimate the SFS for entire chromosome 22 ??
	-> ./realSFS afile.saf.idx -r chr22 

	-> 2) Estimate the 2d-SFS for entire chromosome 22 ??
	-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 

	-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
	-> ./realSFS afile.saf.idx -nSites 500000000 

	-> 4) Estimate the SFS around a gene ??
	-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 

	-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]

	-> See realSFS print for possible print options
	-> Use realSFS print_header for printing the header

	->------------------
	-> NB: Output is now counts of sites instead of log probs!!
	-> NB: You can print data with ./realSFS print afile.saf.idx !!
	-> NB: Higher order SFS can be estimated by simply supplying multiple .saf.idx files!!
	-> NB: Program uses accelerated EM, to use standard EM supply -m 0
```

This command will estimate the SFS for each population:
```
for POP in LWK TSI PEL
do
        echo $POP
        realSFS $POP.saf.idx -P 4 2> /dev/null > $POP.sfs
done
```
The output will be saved in *.sfs files.

You can now have a look at the output file, for instance for the African (LWK) samples:
```
cat LWK.sfs
```
The first value represent the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.

For 20 individuals, how many values do you expect?
```
awk -F' ' '{print NF; exit}' LWK.sfs 
```
Indeed this represents the unfolded spectrum, so it has 2N+1 values with N diploid individuals.
Please note that this maximum likelihood estimation of the SFS should be performed at the whole-genome level to have enough information for the algorithm to converge.
However, for practical reasons, here we could not use large genomic regions (the BAM files we are using only contain reads from a section of chromosome 11).

You can plot the SFS for each pop using this simple R script (included in the ngsTools github website):
```
Rscript $SCRIPTS/plotSFS.R LWK.sfs TSI.sfs PEL.sfs
evince LWK_TSI_PEL.pdf
```

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences.
This can be achieved in ANGSD using:
```
realSFS PEL.saf.idx -bootstrap 10  2> /dev/null > PEL.boots.sfs
cat PEL.boots.sfs
```
This command may take some time.
The output file has one line for each bootstrapped replicate.

It is very useful to estimate a multi-dimensional SFS, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).

An important issue when doing this is to be sure that we are comparing the exact same sites between populations. ANGSD does that automatically and considers only the set of overlapping sites in each population pair. For example, the 2D-SFS for the TSI-PEL population pair is computed as follows:

```
realSFS -P 4 TSI.saf.idx PEL.saf.idx 2> /dev/null > TSI.PEL.sfs
```

The output file is a flattened matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.
```
less -S TSI.PEL.sfs
```
You can plot it, but you need to define how many samples you have per population.
```
Rscript $SCRIPTS/plot2DSFS.R TSI.PEL.sfs 20 20 TSI PEL
evince TSI.PEL.sfs.pdf
```

You can even estimate the SFS with more than 2 populations, for example a 3-D SFS (Do not run this right now, it will take a long time to finish):

```
realSFS -P 4 LWK.saf.idx TSI.saf.idx PEL.saf.idx 2> /dev/null > LWK.TSI.PEL.sfs
```



## Nucleotide diversity

You may be interested in assessing levels of nucleotide diversity within a particular population. We can achieve this using ANGSD by estimating levels of diversity without relying on called genotypes. The SFS is used as a prior to compute allele frequencies probabilities. From these quantities, expectations of various diversity indexes are computed. This can be achieved using the following pipeline, assuming we are computing such indexes for all populations (but separately).

First we compute the allele frequency posterior probabilities and associated statistics (-doThetas) using the SFS as prior information (-pest)
```
for POP in LWK TSI PEL
do
	echo $POP
	angsd -P 4 -b $DATAFOL/$POP.bamlist -ref $REF -anc $ANC -out $POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 \
		-doThetas 1 -pest $POP.sfs &> /dev/null
done
```

Then we are going to index these files and perform a sliding windows analysis using a window length of 50kbp and a step size of 10kbp.
```
for POP in LWK TSI PEL
do
	echo $POP
	# index files
	thetaStat make_bed $POP.thetas.gz &> /dev/null
	# perform a sliding-window analysis
	thetaStat do_stat $POP.thetas.gz -nChr 1 -win 50000 -step 10000 -outnames $POP.thetas &> /dev/null
done
```
Values in this output file are the sum of the per-site estimates for the whole window.

For instance:
```
less -S PEL.thetas.pestPG
```

Finally, you may also be interested in estimating allele frequencies for single SNPs of interest.
In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
Assume that these are the SNPs we are interested in (chromosome and genomic position 1-based):
- 11 61627960 <br>
- 11 61631510 <br>
- 11 61632310 <br>
- 11 61641717 <br>
- 11 61624414 <br>
- 11 61597212 <br>

The file with these positions need to be formatted as (chromosome positions).
```
> $DATAFOL/Data/snps.txt
echo 11 61627960 >> Data/snps.txt
echo 11 61631510 >> Data/snps.txt
echo 11 61632310 >> Data/snps.txt
echo 11 61641717 >> Data/snps.txt
echo 11 61624414 >> Data/snps.txt
echo 11 61597212 >> Data/snps.txt
```
We need to index this file in order for ANGSD to process it. This has already been done in our Data folder using the following command:
```
angsd sites index $DATAFOL/Data/snps.txt
```

We are interested in calculating the derived allele frequencies, so we are using the ancestral sequence to polarize the alleles with '-doMajorMinor 5'.
Note that here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in LWK TSI PEL
do
        echo $POP
        angsd -P 4 -b $DATAFOL/$POP.bamlist -ref $REF -anc $ANC -out $POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites $DATAFOL/Data/snps.txt &> /dev/null
done
```
Inspect the results.
```
zcat LWK.mafs.gz TSI.mafs.gz PEL.mafs.gz
```



## D-statistic (ABBA-BABA test)

You may also be interested in using D-statistics to detect admixture genome-wide. We’ll use 10 individual bam files and look at all possible triplet combinations of the D-statistic. The data has already been downloaded in the $DATAFOL/bams folder, and the bam files have been indexed using samtools:


The list of all the individuals we want to include in our analysis has been prepared as follows:
```
ls $DATAFOL/Data/bams/*.bam > $DATAFOL/abbababatest.bamlist
```

First, we run ANGSD with the option -doAbbababa. ANGSD will use the ancestral file provided via -anc as the outgroup (O). Then, it will calculate D(H1,H2,H3,O) for all possible combinations of H1, H2 and H3 from among the individuals listed in the *bamlist file.
```
angsd -out out -doAbbababa 1 -bam $DATAFOL/abbababatest.bamlist -doCounts 1 -anc $DATAFOL/Data/anc.fa.gz
```
 
Then, we perform a block jack-knife on the results, using a custom R script (this script can be downloaded from the angsd github website).
```
Rscript $SCRIPTS/jackKnife.R file=out.abbababa indNames=/ricco/data/fernando/TutorialFiles/abbababatest.bamlist outfile=D_output
```

The results can be found in the D_output.txt file.
