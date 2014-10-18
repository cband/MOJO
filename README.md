# MOJO 

Minimum Overlap Junction Optimizer (MOJO) is an algorithm developed to identify gene fusions from paired-end transcriptome sequencing data.  The fundamental approach is as follows.  First, clusters of discordant reads are identified by mapping reads to the transcriptome in iterative steps to maximize sensitivity.  Next, candidate fusion junctions are constructed from the exons predicted to be involved in fusions between the pairs of genes.  Reads that cannot be aligned to the canonical transcriptome are mapped to these junctions.  Finally, high confidence fusions are nominated following rigorous filtering steps designed to capture both technical and biological noise.  __Note__: Currently, MOJO can only detect fusions at canonical exon-exon junctions.  

MOJO can be used to detect both somatic fusions from cancer transcriptomes and germline polymorphic fusions in normal tissues.  Reference indices are currently available for human, mouse and drosophila melanogaster.


## Quick Start

There are three main steps to get MOJO up and running:  

1. Download pre-built binaries from the github repository
2. Download one of the pre-built custom reference indices
3. Modify the paths in the template configuration file provided

### Download pre-built binaries

Pre-compiled MOJO binaries for Linux are available for download from [here](http://dmel.uchicago.edu/~chai/MOJO/releases/).  Binaries were generated on Ubuntu (kernel-2.6.32) with GCC 4.8.1 and should be compatible with most 64-bit Linux distributions. 

    > ### Find latest release here: http://dmel.uchicago.edu/~chai/MOJO/releases/
    > wget http://dmel.uchicago.edu/~chai/MOJO/releases/MOJO.<latest_release>.tar.gz
    > tar -zxf MOJO.<latest_release>.tar.gz

    
### Setup environment paths

    > export LD_LIBRARY_PATH=<path_to_mojo_directory>/lib\:$LD_LIBRARY_PATH
	> export PATH=<path_to_mojo_directory>/bin\:$PATH
    
### Download MOJO reference files

Pre-built MOJO custom reference index files are generated for each genome and gene model.  Ensembl/knownGene annotations are transformed into a custom format upon which the spliced/unspliced transcriptome and genome indices (bwa and bowtie2) are built.  This also includes a megablast index of all pairwise comparisons of genes in the corresponding annotation.  Currently, reference indices for seven genome/transcript annotation models are available.  

| Genome      | Transcript model  |   Genes  |   Isoforms  |   Exons  |
| :---------: |:-----------------:| :-------:| :----------:| :-------:|
| hg19        | [TCGA GAF 3.0](http://dmel.uchicago.edu/~chai/MOJO/references/reference.hg19.GAF3.0.tar.gz) **   | 26,627   | 73,900      |  272,808 |
| hg19       | [UCSC knownGene](http://dmel.uchicago.edu/~chai/MOJO/references/reference.hg19.knownGene.tar.gz)    | 30,522 |  78,826  |  279,989   |
| hg19        |[Ensembl](http://dmel.uchicago.edu/~chai/MOJO/references/reference.hg19.Ensembl.tar.gz)           | 57,773 | 196,354  |  568,095   |
| hg38        |[UCSC knownGene](http://dmel.uchicago.edu/~chai/MOJO/references/reference.hg38.knownGene.tar.gz)           | 44,037 | 92,716  |  295,885   |
| mm10     | [UCSC knownGene](http://dmel.uchicago.edu/~chai/MOJO/references/reference.mm10.knownGene.tar.gz)    | 32,182 |  61,396  |  255,134   |
| mm10     | [Ensembl](http://dmel.uchicago.edu/~chai/MOJO/references/reference.mm10.Ensembl.tar.gz)           | 38,924 |  94,545  |  349,676   |
| dm3        | [Ensembl](http://dmel.uchicago.edu/~chai/MOJO/references/reference.dm3.Ensembl.tar.gz)    | 15,681 | 29,172 | 77,793 |


\** ___TCGA GAF 3.0___: General Annotation Format (GAF) is the genomic annotation used by the various analysis working groups of The Cancer Genome Atlas (TCGA) project.  GAF 3.0 is a slightly curated version of UCSC knownGene model. 


    > wget http://dmel.uchicago.edu/~chai/MOJO/reference.<GENEMODEL>.tar.gz
    > tar -zxf reference.<GENEMODEL>.tar.gz -C /mojo_directory/
    
### Setup configuration file for MOJO run

A template configuration file `Sample.configfile.txt` is provided in the top-level directory of MOJO.  The following three parameters need to be configured correctly (all other parameters in the config file are required but can remain unchanged):

    mojo_install_dir       =  <mojo_directory>/bin/
    mojo_reference_dir     =  <mojo_directory>/references.<GENEMODEL>/
    mojo_tools_dir         =  <mojo_directory>/external/


### Running MOJO

    MOJO --config <CONFIG> 
         --sample_name <NAME> 
         --output_dir <OUT_DIR> 
         --fq1 <lane1_1.fastq,lane2_1.fastq,...> 
         --fq2 <lane1_2.fastq,lane2_2.fastq,...> 


### Example

### Test run using reads from K562 cell line ###
[Download](http://dmel.uchicago.edu/~chai/MOJO/test_data/Edgren_KPL4.tar) a test transcriptome comprised of 5 million reads sub-sampled from a K562 transcriptome (PubmedID: 20179022).  
    
    > wget http://dmel.uchicago.edu/~chai/MOJO/test_data/Edgren_KPL4.tar
    > tar -zxvf Edgren_KPL4.tar
    > 
    > MOJO --config hg19.Ensembl.configfile.txt 
           --sample_name Edgren_KPL4 
           --output_dir ./
           --fq1 Edgren_KPL4_1.fastq.gz
           --fq2 Edgren_KPL4_2.fastq.gz
    
###Expected output###

## Installation from source

### Requirements

- A C++11 supported compiler GCC 4.7 or higher is required. Download from: https://gcc.gnu.org/mirrors.html
- CMake 2.8 or higher is required to build MOJO and its dependencies.  The CMake build system can automatically fetch, compile and install dependencies. Download from: http://www.cmake.org/cmake/resources/software.html


### Installing MOJO

MOJO installation requires additional dependencies: `boost-1.55.0`, `bowtie2-2.2.3`, `bwa-0.7.10`, `samtools-1.0.0`, `bamtools-2.3.0`.  A `cmake` makefile is provided to automatically download, compile and install all dependencies.  MOJO has been extensively tested with these specific versions , therefore, it is recommended to presist with these when possible. 

___Note:___ If `boost-1.55.0` is not already installed, this process can take several hours. `cmake` checks the `BOOST_ROOT` environment variable to find boost installation. If you have a local installation of boost-1.55.0, please ensure that this variable is correctly set. 

Download latest MOJO source distribution for the latest release (https://github.com/cband/MOJO-P/releases). Or, just clone the repo with git.

    > tar -zxvf MOJO.<LATEST>-source.tar.gz
    > cd MOJO.<LATEST>-source 

Inside the MOJO source directory, create build directory and run ```cmake``` in it.  

    > mkdir build
    > cd build
    > cmake ..

Run `make`. _(if boost installation is required, this can take a few hours)_

    > make
    > make install

If the build and installation processes are successful, add `libs` and `bin` paths to environment variables LD\_LIBRARY\_PATH and PATH, respectively.  (see __Setup environment paths__ section  above)


Finally, run the following to make sure the installation is successful. See __Example__ above to start a test run.

    > MOJO
    
## MOJO General Usage

### Input parameters


| parameter        | description|
| ------------: | :----------------|
|   --config, -c |  MOJO configuration file  ___[required]___ <br>Template config files (`MOJO.<GeneModel>.configfile.txt`) are provided for each reference. |
|   --output_dir, -o | output directory for the run  ___[required]___  |
|   --sample_name, -s|  sample name  ___[required]___ |
|   --fq1, -1|  comma separated full paths of end 1 lanes  ___[required]___   |
|   --fq2, -2|  comma separated full paths of end 2 lanes   ___[required]___  |
| --cores | number of cores to use [default:  total cpu count]  | 
| --mem, -m | max available memory (in GB) [default: 85% of system memory]. MOJO uses this to appropriately parallelize memory intensive tasks.  | 
| --min_span | threshold for minimum number of discordant reads.  Configured as a function of library size: <br>`ceil(R + X * max(0,ln(LibrarySize/Y))` <br> `R`, `X` and `Y` are constants representing minimum of span reads, a coefficient and a scaling factor (that accounts for library size), respectively.  <br>[default (`R,X,Y`): `--min_span 2,2,80000000`] <br> To set hard cut-offs, set X and Y to 0. |
| --read_through |  a fusion is designated as read-through if the genes are on the same strand and the 5' gene is upstream of the 3' gene, and, the distance between both is less than `read_through`.   [default: 200000] |
|  --junct\_mismatch| max mismatch rate in split-reads aligning to function junctions [default: 0.03]  |


### Output description

Two primary output files are generated to describe the fusion results.  `<sample_name>.fusions` contains a listing of all the fusion calls.  `<sample_name>.fusions.pileup` contains a pileup of reads mapping to all the fusion junctions nominated in `<sample_name>.fusions`.

####<sample_name>.fusions file definition

| ID | column name    | description      |
|:--:| :------------: | :----------------|
| 1. | GeneA_GeneB  | Fusion gene name | 
| 2. | n_discord_AB  | # of discordant read pairs between genes A and B | 
| 3. | n_unique_discord_AB  | # of unique discordant read pairs between A and B  | 
| 4. |  n_anchor_reads | # of anchor reads supporting the fusion junction.  An anchor read is a paired-end read with one end mapping to the fusion junction and the other mapping to either A or B  | 
| 5. | n_high_conf_anchor_reads  | # of high confidence anchor reads.  See __FAQ__ for criteria for a high confidence anchor read | 
| 6. | gene_5p  | 5' gene name | 
| 7. | chrom_5p  | 5' gene chromosome  | 
| 8. | strand_5p  | 5' gene strand | 
| 9. | exon_id_5p  | 5' exon id (identification in MOJO) | 
| 10. | breakpoint_5p  | 5' breakpoint junction position | 
| 11. | breakpoint_region_5p | 5' breakpoint region (CDS, 5'-UTR, 3'-UTR or non-coding RNA) | 
| 12. | gene_3p  | 3' gene name | 
| 13. | chrom_3p  | 3' gene chromosome | 
| 14. | strand_3p  | 3' gene strand | 
| 15. | exon_id_3p  | 3' exon id (identification in MOJO)  | 
| 16. | breakpoint_3p  | 3' breakpoint junction position | 
| 17. | breakpoint_region_3p | 3' breakpoint region (CDS, 5'-UTR, 3'-UTR or non-coding RNA)  | 
| 18. | is_inframe | 0/1 to indicate if the fusion can generate an in-frame transcript | 
| 19. | n_anchor_reads_10bp | # of anchor reads with anchor length betweeen 10-14bp | 
| 20. | n_anchor_reads_15bp | # of anchor reads with anchor length betweeen 15-19bp | 
| 21. | n_anchor_reads_20bp | # of anchor reads with anchor length 20bp or higher | 
| 22. | n_high_conf_anchor_A | # of high confidence anchor reads with the non-junction end mapping to A  | 
| 23. | n_high_conf_anchor_B | # of high confidence anchor reads with the non-junction end mapping to B  | 
| 24. | distance_AB | distance between genes A and B. -1 if inter-chromosomal. 
| 25. | entropy_5p | Entropy of 20bp of junciton sequence from gene A | 
| 26. | entropy_3p | Entropy of 20bp of junction sequence from gene B  | 
| 27. | rpkm_A | RPKM for gene A |
| 28. | rpkm_B | RPKM for gene B |
| 29. | rpkm_A_5p | RPKM for 5' fragment of gene A |
| 30. | rpkm_A_3p | RPKM for 3' fragment of gene A |
| 31. | rpkm_B_5p | RPKM for 5' fragment of gene B |
| 32. | rpkm_B_3p | RPKM for 3' fragment of gene B |
| 33. | n_concords_A_5p | # of concordant reads mapping to region of gene A that is 5' of the fusion junction (_region involved in fusion_)| 
| 34. | n_concords_A_3p | # of concordant reads mapping to region of gene A that is 3' of the fusion junction | 
| 35. | n_concords_B_5p | # of concordant reads mapping to region of gene B that is 5' of the fusion junction  | 
| 36. | n_concords_B_3p | # of concordant reads mapping to region of gene B that is 3' of the fusion junction (_region involved in fusion_) | 
| 37. | n_concords_span_AA | # of concordant reads with the ends spanning the breakpoint in gene A | 
| 38. | n_concords_span_BB | # of concordant reads with the ends spanning the breakpoint in gene B | 
| 39. | n_concords_junct_AA | # of concordant reads mapping to gene A with one end mapping to the breakpoint (canonical exon-exon junction of gene A)  | 
| 40. | n_concords_junct_BB | # of concordant reads mapping to gene B with one end mapping to the breakpoint (canonical exon-exon junction of gene B) | 
| 41. | n_discords_A5p_B3p | # of discordant reads with one end mapping to region that is 5' of the breakpoint in gene A and the other end mapping to region that is 3' of the breakpoint in gene B | 
| 42. | n_discords_A3p_B5p | # of discordant reads with one end mapping to region that is 3' of the breakpoint in gene A and the other end mapping to region that is 5' of the breakpoint in gene B | 
| 43. | n_discords_A5p_B5p | # of discordant reads with one end mapping to region that is 5' of the breakpoint in gene A and the other end mapping to region that is 5' of the breakpoint in gene B | 
| 44. | n_discords_A3p_B3p | # of discordant reads with one end mapping to region that is 3' of the breakpoint in gene A and the other end mapping to region that is 3' of the breakpoint in gene B | 
| 45. | coding_sequences | predicted fusion ORFs.  Format: `[isoformID_A-isoformID_B]:sequence`|
| 46. | transcribed_sequences | potential transcribed fusion transcripts. `[isoformID_A-isoformID_B]:sequence`  | |

####<sample_name>.fusions.pileup file definition

A pileup of anchor reads mapping to the fusion junction.  

## FAQs

### Installation issues  ###

__1. `xxxx.h: No such file or directory` errors during compiling from source__

This error can be due to a missing dependency that is expected to be already available on most Unix distributions.  Please ensure that the following are installed before posting an issue on github.

    sudo apt-get install libpng-dev         ## required by Blat
    sudo apt-get install python-dev         ## required by Boost
    sudo apt-get install libbz2-dev         ## required by Boost
    sudo apt-get install libncurses5-dev    ## required by Samtools


__2. `cannot find -lcurses` error during samtools installation__

In the `samtools` makefile `./external/source/samtools/Makefile`, change 
     
     LIBCURSES= -lcurses

to

     LIBCURSES= -lncurses
     

See: http://seqanswers.com/forums/showthread.php?t=6669

### General Usage ### 

__0. Some definitions__

An ___anchor read___ is a paired-end read with one end mapping to the fusion junction (___split-end read___) and the other end (___other-end read___) mapping to one of the two genes of the fusion pair.  And, the ___anchor length___ is the minimum overhang of the split-end read.  Anchor length is one of the key determinants of the specificity of a fusion call.   

__1. Why is anchor length not provided as a parameter to configure?__

MOJO implicitly uses a minimum anchor length threshold of 10bps and reports all anchor reads with lengths >= 10bps.  Numbers of anchor reads with length 10, 15 and 20bps are reported in results allowing for customized filtering. For more precise information on all the anchor lengths, see .fusions.pileup file.

__2. What is a 'high confidence' anchor read (column number 5 of output)?__

An anchor read is classified as high confidence if it satisfies all of the following criteria: 

1) The split-end read does not have alternate alignments to the genome or transcriptome, 

2) The other-end read maps uniquely to one of the two genes in the fusion pair, 

3) The split-end read has anchor-length >= 20bps, or, the gene to which the smaller overhang of the split-read maps to is also the gene to which the other-end maps to.   

__3. How to interpret the confidence level in a fusion call?__

The number of `n_high_conf_anchor_reads` is strongly associated with the confidence level in a fusion call.  We observe that `n_high_conf_anchor_reads >= 2` will yield high specificity.  If sequencing depth is low or low expressed fusions are of interest, then consider lowering this threshold.

In addition, the columns `n_high_conf_anchor_A` and `n_high_conf_anchor_B` contain the number of high confidence anchor-reads with the other-ends mapping to gene A and gene B, respectively.  Both columns being non-zero is a strong indication for a high confidence call.  However, if sequencing depth is limited or if the breakpoint is at the ends of transcripts (due to 5'/3' coverage depletion in RNA-seq), either of the two columns can be 0.  


__4. How does MOJO perform in comparison with other published methods?__

Comparisons of MOJO with other published methods will soon be available on the Wiki. 


