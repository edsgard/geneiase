
###About
GeneiASE is a software for detection of condition-dependent allele
specific expression in single individuals. GeneiASE does not require
haplotype phasing and performs consistently over a range of read
depths and ASE effect sizes. See the <a href="http://www.nature.com/articles/srep21134">paper</a> for further
information.

Copyright © 2015 Daniel Edsgärd, Olof Emanuelsson<br>
GeneiASE is available free to use, under the <a href="./LICENSE">GNU GPL version 3 license</a>.

###Download
You find the latest stable release <a href="https://github.com/edsgard/geneiase/releases">here</a>.

###Installation
#####Prerequisites
There are a number of R package dependencies. The dependencies can be installed from within R by:
 ```R
 install.packages(c('getopt', 'binom', 'VGAM'))
```

#####Test the installation
GeneiASE can be run from the shell prompt by entering the downloaded
and unzipped directory (geneiase) and then executing the program, residing in the bin directory:
 ```Shell
cd geneiase
bin/geneiase -t static -i test/static.test.input.tab -b 100
bin/geneiase -t icd -i test/icd.test.input.tab -b 100
```

#####Optionally add geneiase to your shell PATH
By adding the 'geneiase/bin' directory to your shell PATH you can execute
the program without needing to entering the directory where it resides or using the full path.


###Running GeneiASE

#####Required arguments
Only two arguments are required.<br>
 1. -t, followed by a string with the allowed values "static" or "icd", specifying if static or individual condition-dependent ASE is to be tested for.
 2. -i, followed by the input file name. The input file should contain tab-separated columns. In the case of static ASE there should be four columns: feautureID, snpID, alternative allele count, reference allele count, and in the case of icd-ASE there should be six columns: feautureID, snpID, Untreated alternative allele count, Untreated reference allele count, Treated alternative allele count, Treated reference allele count.
 
See **Test the installation** above for an example.

##### Description of arguments
For detailed help on each available argument run the program with the -h flag:
```Shell
geneiase -h
```

##### Input files
Examples of how to format the input files you find in the directory 'test' as part of the code-bundle.

#####Output description
 A test-statistica, s, is generated for each variant reflecting the effect-size (see <a href="http://www.nature.com/articles/srep21134">paper</a>). In a meta-analysis approach the effect sizes for all variants within a gene is combined. We provide several simple gene-wise measures based on the variant effect-sizes. The p-value is based on the Liptak-Stouffer method (column liptak.s) and the generation of a null distribution from resampling of a parametric model (see <a href="http://www.nature.com/articles/srep21134">paper</a>). The fields in the tab-separated output file are:
- feat: FeatureID as specified in the input file (typically a gene identifier)
- n.vars: Number of variants within the gene
- mean.s: Mean of s across the variants within the gene
- median.s: Median of s across the variants within the gene
- sd.s: Standard deviation of s across the variants within the gene
- cv.s: Coefficient of variation of s across the variants within the gene
- liptak.s: Stouffer-Liptak combination of s
- p.nom: Nominal p-value
- fdr: Benjamini-Hochberg corrected p-value

###Typical workflow for generation of the input table containing allelic read counts
 We are not aware of any pipeline that on a genome-wide scale is designed for the particular purpose of generating a table with allelic read counts. The generation of allelic counts starting from raw reads can be specific to sequencing technology and desired tools (such as choice of read mapping and variant calling software). However, to facilitate the generation of such a table we list a possible workflow below. 
 1. Quality control of reads, such as trimming (FASTQ)
 2. Map reads (BAM)
 3. Quality control mapped reads, such as PCR duplicate removal
 4. Call variants (VCF)
 5. Filter variants: Heterozygosity, read depth, and satisfying other mapping quality and variant calling criteria
   5.2. Optional: Assess mapping bias for each variant and either use the mapping bias estimate of each variant in downstream analysis or filter out variants exhibiting mapping bias.
   5.3. Optional: If there are several samples from the same individual, such as in the case of cd-ASE, filter out variants which differ between the samples.
 6. Annotate variants
 7. Filter on relevant annotation, such as within a gene or presence in dbSNP.
 8. Get allelic counts for filtered variants based on a pileup of mapped reads.

### Citing GeneiASE
If you use GeneiASE, please cite it as follows:
Edsgärd D. et al., GeneiASE: Detection of condition-dependent and static allele-specific expression from RNA-seq data without haplotype information, <em>Scientific Reports</em>, 2016

### Contact
Please contact the corresponding author if you have technical issues or other comments or questions:<br/>
<a href="mailto:olofem (at) kth (dot) se">Olof Emanuelsson</a>


