Exercises on demographic analyses with ngsTools and ANGSD - Copenhagen 2019 course
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

We will use the same data that we used yesterday during Matteo's tutorial. As a reminder, this includes 50 BAM files of a section of chr2 from human samples (i.e. 10 African, 10 Native American, 10 European, 10 East Asian, and 10 Latino individuals), a reference genome, and putative ancestral sequence.

The data can be found here:
/ricco/data/matteo/Data

As before, let’s create a shortcut for this directory, as we will refer to it quite a lot in our command lines:
```
DATA=/ricco/data/matteo/Data
```

We also have a reference and an ancestral sequence in FASTA format. Let’s also make shortcuts for them:
```
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz
```

Side-note: in case an ancestral sequence is not available, analyses on the SFS and nucleotide diversity can be carried out using the reference sequence to polarize your data. Please be aware that, under this scenario, some quantities (e.g. the unfolded joint site frequency spectrum) will be nonsense. Please also note that, since we are randomly subsampling reads here, your results in this tutorial may (slightly) differ from what is written here. 

Throughout the tutorial, we will also use a number of R scripts that can be downloaded from the ngsTools github website. Let’s make a shortcut to them too:
```
SCRIPTS=/ricco/data/fernando/Scripts
```

Summary statistics using ANGSD
-----------------------------------

Nucleotide diversity indexes and measures of population differentiation can be estimated taking data uncertainty into account both with ANGSD and ngsTools.

## The Site Frequency Spectrum (SFS)

One of the most important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS). 
SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. 
The SFS is informative about the demography of the population or about selective events (when calculated at a local scale).

We will use ANGSD to estimate the SFS with an example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
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
We need to slightly modify the filtering options as each population has 10 samples. 
So now we set `-minInd 10 -setMinDepth 20 -setMaxDepth 200`.
Also, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarize our alleles (to ancestral and derived states).

We cycle across all populations:
```
for POP in AFR EUR LAT EAS NAM
do
	echo $POP
	angsd -P 4 -b $DATA/$POP.bams -ref $REF -anc $ANC -out $POP \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
		-GL 1 -doSaf 1 &> /dev/null
done
```

Some basic filtering consists in removing, for instance, reads with low quality and/or with multiple hits, and this can be achieved using the parameters ```-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1```. As input we give the list of BAM files with option `-b` and then specify the references sequence with `-ref` and the prefix for output files with `-out`. Additionally, ```-C 50``` reduces the effect of reads with excessive mismatches, while ```-baq 1``` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) to rule out false SNPs close to INDELS, and ```-trim 0``` means that we are not trimming the ends of reads. With ```-minMapQ 20``` we filter out reads with low mapping quality. Finally, ```-P 4``` means that I am using 4 threads.

Have a look at the output file.
```
realSFS print LAT.saf.idx | less -S
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
for POP in AFR EUR LAT EAS NAM
do
        echo $POP
        realSFS $POP.saf.idx -P 4 2> /dev/null > $POP.sfs
done
```
The output will be saved in *.sfs files.

You can now have a look at the output file, for instance for the African (AFR) samples:
```
cat AFR.sfs
```
The first value represent the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.

For 10 diploid individuals, how many values do you expect?
```
awk -F' ' '{print NF; exit}' AFR.sfs 
```
Indeed this represents the unfolded spectrum, so it has 2N+1 values with N diploid individuals.
Please note that this maximum likelihood estimation of the SFS should be performed at the whole-genome level to have enough information for the algorithm to converge.
However, for practical reasons, here we could not use large genomic regions (the BAM files we are using only contain reads from a section of chromosome 2).

You can plot the SFS for each pop using this simple R script (included in the ngsTools github website):
```
Rscript $SCRIPTS/plotSFS.R AFR.sfs EUR.sfs LAT.sfs
evince AFR_EUR_LAT.pdf
```

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences.
This can be achieved in ANGSD using:
```
realSFS LAT.saf.idx -bootstrap 10  2> /dev/null > LAT.boots.sfs
cat LAT.boots.sfs
```
The output file has one line for each bootstrapped replicate.

It is very useful to estimate a multi-dimensional SFS, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).

An important issue when doing this is to be sure that we are comparing the exact same sites between populations. ANGSD does that automatically and considers only the set of overlapping sites in each population pair. For example, the 2D-SFS for the EUR-AFR population pair is computed as follows:

```
realSFS -P 4 EUR.saf.idx AFR.saf.idx 2> /dev/null > EUR.AFR.sfs
```

The output file is a flattened matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.
```
less -S EUR.AFR.sfs
```
You can plot it, but you need to define how many samples you have per population. Again, we use a handy script that can also be obtained from the ngsTools github website:
```
Rscript $SCRIPTS/plot2DSFS.R EUR.AFR.sfs 10 10 EUR AFR
evince EUR.AFR.sfs.pdf
```

You can even estimate the SFS with more than 2 populations, for example a 3-D SFS:

```
realSFS -P 4 AFR.saf.idx EUR.saf.idx LAT.saf.idx 2> /dev/null > AFR.EUR.LAT.sfs
```

Side-note: when performing demographic inference, estimating the SFS from the sample is one of the first steps in the analyses. Afterwards, one may want to run downstream programs (like dadi or fastsimcoal2) that attempt to infer the best demographic model that can be fitted to the SFS.

## Nucleotide diversity

You may be interested in assessing levels of nucleotide diversity within a particular population. We can achieve this using ANGSD by estimating levels of diversity without relying on called genotypes. The SFS is used as a prior to compute allele frequencies probabilities. From these quantities, expectations of various diversity indexes are computed. This can be achieved using the following pipeline, assuming we are computing such indexes for all populations (but separately).

First we compute the allele frequency posterior probabilities and associated statistics (-doThetas) using the SFS as prior information (-pest)
```
for POP in AFR EUR LAT EAS NAM
do
	echo $POP
	angsd -P 4 -b $DATA/$POP.bams -ref $REF -anc $ANC -out $POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 -doThetas 1 -pest $POP.sfs &> /dev/null
done
```

Then we are going to index these files and perform a sliding windows analysis using a window length of 50kbp and a step size of 10kbp.
```
for POP in AFR EUR LAT EAS NAM
do
	echo $POP
	# perform a sliding-window analysis
	thetaStat do_stat $POP.thetas.idx -win 50000 -step 10000 -outnames $POP.thetas
done
```
Values in this output file are the sum of the per-site estimates for the whole window.

For instance:
```
less -S LAT.thetas.pestPG
```

The .pestPG file is a 14 column file (tab-separated). The first column contains information about the region. The second and third column are the reference name and the center of the window. We then have 5 different estimators of theta, these are: Watterson's estimator, Tajima's estimator (pi), Fu&Li's estimator, Fay's H and L. And then we have 5 different neutrality test statistics: Tajima's D, Fu&Li's F, Fu&Li's D, Fay's H and Zeng's E. The final column is the effetive number of sites with data in the window.

## Estimating per-SNP allele frequencies from low-coverage data

Using ANGSD, we can also estimate allele frequencies for specific SNPs of interest. In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
Assume that these are the SNPs we are interested in (chromosome and genomic position 1-based):
- 2 109000112 <br>
- 2 109000319 <br>
- 2 109000331 <br>
- 2 109000433 <br>
- 2 109000373 <br>
- 2 109000433 <br>

The file with these positions need to be formatted as (chromosome positions). To create the file, we run the following commands.

```
echo 2 109000112 > snps.txt
echo 2 109000319 >> snps.txt
echo 2 109000331 >> snps.txt
echo 2 109000433 >> snps.txt
echo 2 109000373 >> snps.txt
echo 2 109000433 >> snps.txt
```

We need to index this file in order for ANGSD to process it:
```
angsd sites index snps.txt
```

We are interested in calculating the derived allele frequencies, so we are using the ancestral sequence to polarize the alleles with '-doMajorMinor 5'.
Note that here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in AFR EUR LAT EAS NAM
do
        echo $POP
        angsd -P 4 -b $DATA/$POP.bams -ref $REF -anc $ANC -out $POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites snps.txt &> /dev/null
done
```
Inspect the results.
```
zcat AFR.mafs.gz EUR.mafs.gz LAT.mafs.gz
```



## D-statistic (ABBA-BABA test)

Finally, you may also be interested in using D-statistics to detect admixture genome-wide. We’ll use 10 individual bam files and look at all possible triplet combinations of the D-statistic. The data has already been downloaded in the /ricco/data/fernando/TutorialFiles/Data/bams folder, and the bam files have been indexed using samtools:

Let's make a shortcut for the folder where the files are located:

```
DFILESFOL=/ricco/data/fernando/TutorialFiles/Data/bams
```

The list of all the individuals we want to include in our analysis should be in a list file:
```
ls $DFILESFOL/*.bam > abbababatest.bams
```

First, we run ANGSD with the option -doAbbababa. ANGSD will use the ancestral file provided via -anc as the outgroup (O). Then, it will calculate D(H1,H2,H3,O) for all possible combinations of H1, H2 and H3 from among the individuals listed in the *bams file.
```
angsd -out out -doAbbababa 1 -bam abbababatest.bams -doCounts 1 -anc $DFILESFOL/anc.fa.gz
```
 
Then, we perform a block jack-knife on the results, using a custom R script (this script can be downloaded from the angsd github website).
```
Rscript $SCRIPTS/jackKnife.R file=out.abbababa indNames=abbababatest.bams outfile=D_output
```

The results can be found in the D_output.txt file:

```
less D_output.txt
```

In this file, H1, H2 and H3 are the 3 individuals in the tree that are not the outgroup. H1 and H2 are the "test" individuals and H3 is the potential introgressor.

nABBA the total counts of ABBA patterns

nBABA the total counts of BABA patterns

Dstat is test statistic: (nABBA-nBABA)/(nABBA+nBABA). A negative value means that H1 is closer to H3 than H2 is. A positive value means that H2 is closer to H3 than H1 is.

The JackEst column is another estimate of the abbababa statistic that is bias-corrected. This value should be extremely similar to the value in the Dstat column.

SE is the estimated m-delete blocked Jackknife standard error of the estimate used to obtain the Z-value.

Z is the Z-value that can be used to determine the significance of the D-statistic test. An absolute value of the Z score above 3 is often used as a critical value. However, note that this does not take into account the fact that we are performing multiple tests.
