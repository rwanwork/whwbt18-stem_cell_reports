Stem Cell Reports 2018 Workflow
===============================

Introduction
------------

This repository, **scr2018**, consists of a set of programs for performing the analysis for this manuscript:

    E. van der Wal, P. Herrero-Hernandez, R. Wan, M Broeders, S. L M In 't Groen, T. J. M. van Gestel, W. F. J. van IJcken, T. H. Cheung, A. T. van der Ploeg, G. J. Schaaf, W. W. M. Pim Pijnappel.  Large-Scale Expansion of Human iPSC-Derived Skeletal Muscle Cells for Disease Modeling and Cell-Based Therapeutic Strategies, 10(6), pg. 1975 - 1990, Stem Cell Reports, 2018.


Overview
--------

This workflow does the following:

1.  Download and rename data that has been processed with Hisat2 and Ballgown (single-pass, no novel transcripts).
2.  Apply Ballgown to create a data table of expressed genes.  Ballgown employs the following filtering criteria:
    * Select rows where the row mean is greater than 1 (i.e., remove genes where expression levels are low throughout all data sets).
    * Select rows where at least one cell value is greater than the threshold.
3.  Create the following figures:
    * Hierarchical clustering with 3 different distance measures
    * Hierarchical clustering heatmap with 3 different distance measures
    * Principal Component Analysis (PCA)
    * t-Distributed Stochastic Neighbor Embedding (t-SNE)
4.  For two given samples, average across their replicates.  It is assumed that each sample is of the form "XXX_D" somewhere in the sample name, where XXX is the sample name and D is a number.
5.  An [MA Plot](https://en.wikipedia.org/wiki/MA_plot) is created.  The "M" is defined as "log_2 (sample1 / sample2)".
6.  For the set of significantly expressed genes (based on adjusted p-value and fold change cut-offs), the up- and down-regulated genes are separated into two sets and sent to [DAVID](https://david.ncifcrf.gov/) for gene enrichment.
7.  Up- and down-regulated genes are plotted as two horizontal bar graphs.  (The left figure is up-regulated; the right figure is down-regulated.)

Notes about the "all" rule:
* If the reciprocal of the calculation of "M" is desired, then in the "all" rule in Snakefile.py, reverse the order of the two samples.  i.e., Change "Results/combo7/MPC-ASC.done" to "Results/combo7/ASC-MPC.done".
* The "all" rule can be modified if two samples do not need to be compared at all.  This is true when only the global heatmap is required.  This will eliminate any accesses to the DAVID pipeline and save at least an hour per combination.  To enable this, add rules such as:  "Results/comboX/ballgown.done" in Snakefile.py, where X is the combination number.

The output of Ballgown will differ based on which samples are included in the analysis.  For example, a Ballgown run with Samples A and B will yield different results for these two samples compared to another Ballgown run with Samples A, B, and C.  Their individual FPKM values will be the same, but, as shown in step #2 of the workflow, different genes will be removed.  Five combinations of samples have been created.  They are:

1.  All Pompe samples (8).
2.  Concatenation of the technical replicates (4).
3.  All Pompe samples plus Tom's 6 samples (14).
4.  Concatenation of the technical replicates with Tom's samples (10).
7.  All samples (excluding the concatenated samples).

It is suggested that the concatenated samples (#2 and #4) be ignored.  Previously (about 6 months ago), combinations 5 and 6 involved random sampling of reads.  They are no longer performed.


Preliminaries
-------------

1.  Clone this repository.
2.  Create an account with [DAVID](https://david.ncifcrf.gov/content.jsp?file=WS.html).
3.  Place the registered e-mail address in the variable "EMAIL" of the file global-vars.py .
4.  Expand the data archive in the main directory so that a directory called Data/ is created.  The directory structure should look something like this:

    .  
    ├── combos  
    ├── Data  
    ├── Perl  
    └── R ── modules  
        
5.  Type the following commands:

    touch Data/combo1.download Data/combo1.rename  
    touch Data/combo2.download Data/combo2.rename  
    touch Data/combo3.download Data/combo3.rename  
    touch Data/combo4.download Data/combo4.rename  
    touch Data/combo7.download Data/combo7.rename  
    
    These files are used to indicate that the data files have been downloaded and renamed.  As all of this has been done for the data.tar.gz archive, we want to skip these two steps (the Download and Rename_Data rules in Snakefile.py ).


Setup
-----

1.  Download the latest version of Anaconda from https://www.anaconda.com/download/#linux .  Select the one that uses Python 3.x .
2.  Install Anaconda according to the web site's instructions.
3.  Create an environment called "pompe" by typing the following command:  conda env create -f environment.yml
4.  Switch to this environment by typing:  source activate pompe
5.  Go the directory where "Snakefile.py" exists and run the pipeline by typing:  snakemake --snakefile Snakefile.py -p .  Due to potential [race conditions](https://en.wikipedia.org/wiki/Race_condition), it is advised that only one core is used.  Additional cores can be attempted with the --cores option; at worse, an error should result.


Parameters
----------

*  Placed the registered e-mail for DAVID in the variable EMAIL of global-vars.py .
*  The genes for the gene heatmap are defined in R/modules/about-data.R in the HEATMAP_GENES variable.


Sample Run
----------

See the attached 20180209-typescript.txt for a transcript of a sample run using a single core.  The run time on a server was 42 minutes.


Sample Run Notes
----------------

In the rule "Average_Replicates", check that this output is correct:

    [1] "Number of samples for MPC : 4"
    [1] "MPC_1A" "MPC_1B" "MPC_2A" "MPC_2B"
    [1] "Number of samples for ASC : 2"
    [1] "ASC_1" "ASC_2"

These are the two experimental conditions and they should be averaged correctly.  In this output, four samples are averaged across for "MPC" while two samples are averaged across for "ASC".

Also, in the rule "Plot_MAPlot", there is the phrase:  "Set of genes which did not have a matching Entrez ID:".  This is a list of genes which were dropped from gene enrichment because they could not be mapped from gene symbol (alphanumeric) to Entrez ID (numerical).

In a single run of the full workflow, two accesses to the DAVID API are made.  A5 minute pause added in so that the web server gets a break as frequent accesses may ban further accesses from your computer.


Mapping of SRR to GSE/GSM
-------------------------

Separate from the Snakemake workflow, is the `Perl/ncbi-sample-table.pl` script which takes the `config.yaml` file and matches GSE and GSM accession numbers to the SRR read identifiers.  In order to do this, download this file [1], as explained in this forum [2].

Then, run the script as follows:  `Perl/ncbi-sample-table.pl <config.yaml >sample-table.tsv`

This data table has a mix of information.  One can see this by looking at column 7:

    $ cut -f 7 SRA_Accessions.tab |
    sort | uniq -c
     190461 ANALYSIS
     5155776 EXPERIMENT
     8235787 RUN
     4072356 SAMPLE
     157042 STUDY
     968818 SUBMISSION
          1 Type

The script performs the following look-up:
  SRR --> SRX (experiment) --> GSM (experiment alias)
  SRR --> SRA (submission) --> GSE (submission alias, but with "GEO: " removed)
          
This takes about 5 minutes to run.

[1]  ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
[2]  https://www.biostars.org/p/244150/


Limitations
-----------

See [DAVID's documentation](https://david.ncifcrf.gov/content.jsp?file=WS.html) for an explanation of the following:

*  DAVID's API only allows at most 3,000 genes.  If more are found in Perl/david-chart-report.pl, then the variable PVALUE_THRESHOLD in R/maplot.R needs to be lowered.
*  DAVID's API only allows at most 200 accesses per day from a single user or computer.  Register additional addresses (i.e., of collaborators, with their permission), to get around the first restriction.  To get around the second one, re-run the script from another computer.  Alternatively, go home and sleep and continue working the next day.  :-)

The filtering criteria used means that the set of genes when comparing ASC and MPC (for example) will be different depending on which combination is used.

Also, the bottleneck of the entire pipeline is the access to DAVID.  The time needed to run the entire pipeline could be an hour for each item in the "all" rule of Snakefile.py .

    
Known Bugs
----------

1.  Type `source activate pompe` before running Snakemake.  Even though this command has been included in every rule, it doesn't seem to work due to some bug with Conda and how it checks for uninitialized variables.  Perhaps this will be fixed in a later version.
2.  In the rule "all" of Snakefile.py, the averaging of replicates is specified by providing two samples, separated by a hyphen.  This is because all samples use an underscore and never a hyphen in their names.  R/avg-replicates.R assumes that every sample name is followed by "_\d" (i.e., an underscore and a digit like "ASC_1").  This works for the samples we are interested in for this manuscript, but not every sample can be described in this way.  For example, the "CMC" and "SMC" samples are exceptions.  Do check that they are correctly aggregated!
3.  Sometimes Perl/david-chart-report.pl will exit with this error message:

    Uncaught exception from user code:
      502 Proxy Error at Perl/david-chart-report.pl line 217.
       SOAP::Lite::__ANON__(SOAP::Lite=HASH(0x25bdf28), "\x{a}syntax error at line 1, column 49, byte 49 at /expt/anaconda"...) called at /expt/anaconda3/envs/pompe/lib/perl5/site_perl/5.22.0/SOAP/Lite.pm line 3832
       SOAP::Lite::call(SOAP::Lite=HASH(0x25bdf28), "authenticate", "XXX") called at /expt/anaconda3/envs/pompe/lib/perl5/site_perl/5.22.0/SOAP/Lite.pm line 3792
       SOAP::Lite::__ANON__(SOAP::Lite=HASH(0x25bdf28), "XXX") called at Perl/david-chart-report.pl line 217
       
    The cause of this problem is unknown but one guess is that there are too many accesses to DAVID within a certain amount of time.  One solution is to modify the "all" rule so that one combination is run at a time.



Further Information
-------------------

See this [page](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) for further information about Conda.


