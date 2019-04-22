#!/usr/bin/perl
#  Author:  Raymond Wan
#  Organizations:  Division of Life Science, 
#                  Hong Kong University of Science and Technology
#  Copyright (C) 2016-2017, Raymond Wan, All rights reserved.

use FindBin qw ($Bin);
use lib $FindBin::Bin;  ##  Search the directory where the script is located

use diagnostics;
use strict;
use warnings;

use AppConfig;
use AppConfig::Getopt;
use Pod::Usage;

use SOAP::Lite;
use HTTP::Cookies;


########################################
##  Important constants
########################################

my $EASE_THRESHOLD = 0.1;
my $GENE_COUNT = 2;


########################################
##  Important variables
########################################

##  Arguments provided by the user
my $verbose_arg = 0;
my $summary_arg = 0;
my $debug_arg = 0;
my $email_arg = "";
my $regulate_arg = "";

my %id_to_name;
my @ids;


########################################
##  Important functions
########################################

sub MapToNames {
  my ($genes) = @_;

  $genes =~ s/ //gs;
  my @new_array;
  my @array = split /,/, $genes;
  for (my $k = 0; $k < scalar (@array); $k++) {
    push (@new_array, $id_to_name{$array[$k]});
  }
  $genes = join (",", @new_array);
  
  return ($genes);
}


########################################
##  Process arguments
########################################
##  Create AppConfig and AppConfig::Getopt objects
my $config = AppConfig -> new ({
           GLOBAL => {
             DEFAULT => undef,     ##  Default value for new variables
           }
  });

my $getopt = AppConfig::Getopt -> new ($config);

$config -> define ("email", {
           ARGCOUNT => AppConfig::ARGCOUNT_ONE,
           ARGS => "=s"
  });                        ##  E-mail address to use
$config -> define ("regulate", {
           ARGCOUNT => AppConfig::ARGCOUNT_ONE,
           ARGS => "=s"
  });                        ##  Up or down regulated genes?
$config -> define ("verbose!", {
           ARGCOUNT => AppConfig::ARGCOUNT_NONE
  });                        ##  Verbose output
$config -> define ("summary!", {
           ARGCOUNT => AppConfig::ARGCOUNT_NONE
  });                        ##  Summary output
$config -> define ("debug!", {
           ARGCOUNT => AppConfig::ARGCOUNT_NONE
  });                        ##  Debug output
$config -> define ("help!", {
           ARGCOUNT => AppConfig::ARGCOUNT_NONE
  });                        ##  Help screen

##  Process the command-line options
$config -> getopt ();


########################################
##  Validate the settings
########################################

if ($config -> get ("help")) {
  pod2usage (-verbose => 0);
  exit (1);
}

$verbose_arg = 0;
if ($config -> get ("verbose")) {
  $verbose_arg = 1;
}

$summary_arg = 0;
if ($config -> get ("summary")) {
  $summary_arg = 1;
}

$debug_arg = 0;
if ($config -> get ("debug")) {
  $debug_arg = 1;
}

if (!defined ($config -> get ("email"))) {
  printf STDERR "EE\tThe option --email requires an e-mail address.\n";
  exit (1);
}
$email_arg = $config -> get ("email");

if (!defined ($config -> get ("regulate"))) {
  printf STDERR "EE\tThe option --regulate requires one of 'up', 'down', or 'combined' for M values greater than [or equal to] 0, less than 0, or combined, respectively.\n";
  exit (1);
}
$regulate_arg = $config -> get ("regulate");

if (($regulate_arg ne "up") && ($regulate_arg ne "down") && ($regulate_arg ne "combined")) {
  printf STDERR "EE\tThe option --regulate requires one of the word 'up', 'down', or 'combined'.\n";
  exit (1);
}


########################################
##  Summarize the settings
########################################

if ($verbose_arg) {
  printf STDERR "II\tMaximum EASE threshold:  %f\n", $EASE_THRESHOLD;
  printf STDERR "II\tMinimum genes for each term:  %u\n", $GENE_COUNT;
  printf STDERR "II\tRegistered e-mail address:  %s\n", $email_arg;
  printf STDERR "II\tRegulate:  %s\n", $regulate_arg;
}


########################################
##  Read in the gene list
########################################

my $skipped_genes = 0;
my $up_regulated = 0;
my $down_regulated = 0;
my $total_genes = 0;

my $header = <STDIN>;
while (<STDIN>) {
  my $line = $_;
  chomp ($line);
  $total_genes++;
  
  my ($gene_name, $M, $gene_id) = split /\t/, $line;
  
  if ($gene_id eq "NA") {
    $skipped_genes++;
    next;
  }

  if ($M < 0) {
    $down_regulated++;
  }
  
  if ($M >= 0) {
    $up_regulated++;
  }
  
  ##  If we want up-regulate and $M is less than 0, skip this gene
  if (($regulate_arg eq "up") && ($M < 0)) {
    next;
  }
  
  ##  If we want down-regulate and $M is greater than or equal to 0, skip this gene
  ##    Thus, "combined" will not be skipped
  if (($regulate_arg eq "down") && ($M >= 0)) {
    next;
  }
  
  $id_to_name{$gene_id} = $gene_name;
  
  push (@ids, $gene_id);
}
my $ids_as_str = join (",", @ids);

if ($verbose_arg) {
  printf STDERR "II\tTotal number of genes read:  %u\n", $total_genes;
  printf STDERR "II\tSkipped [unmapped] genes:  %u\n", $skipped_genes;
  printf STDERR "II\tUp-regulated genes:  %u\n", $up_regulated;
  printf STDERR "II\tDown-regulated genes:  %u\n", $down_regulated;
  printf STDERR "II\tNumber of mapped genes:  %u\n", scalar (@ids);
  printf STDERR "II\tGenes:  %s\n", $ids_as_str;
}


########################################
##  Authenticate
########################################

my $soap = SOAP::Lite                             
  -> uri ('http://service.session.sample')                
  -> proxy ('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService', cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

my $check = $soap -> authenticate ($email_arg) -> result;
if (lc ($check) ne "true") {
  printf STDERR "Could not authenticate user %s!\n", $email_arg;
  exit (1);
}


########################################
##  Miscellaneous functions to gather information
########################################

# my $conversionTypes = $soap -> getConversionTypes () -> result;
# printf STDERR "\nConversion Types: \n$conversionTypes\n"; 
   
# my $allCategoryNames= $soap ->getAllAnnotationCategoryNames()->result;         
# printf STDERR "\nAll available annotation category names:  %s\n", $allCategoryNames;
 

########################################
##  Add a new list
########################################

my $inputIds = $ids_as_str;
my $idType = 'ENTREZ_GENE_ID';
my $listName = 'Genes';
my $listType = 0;  #  to add background list, set listType=1
my $list = $soap -> addList ($inputIds, $idType, $listName, $listType) -> result;
if ($verbose_arg) {
  printf STDERR "II\t%s of list was mapped\n", $list;
}

#list all species  names
my $allSpecies= $soap -> getSpecies () -> result;
if ($verbose_arg) {
  printf STDERR "II\tAll species:  %s\n", $allSpecies;
}

#list current species  names
my $currentSpecies= $soap -> getCurrentSpecies () -> result;         
if ($verbose_arg) {
  printf STDERR "II\tCurrent species:  %s\n", $currentSpecies;
}


########################################
##  Set categories
########################################

#set user defined categories 
#my $categories = $soap ->setCategories("BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,UP_SEQ_FEATURE")->result;
#to user DAVID default categories, send empty string to setCategories():
# my $categories = $soap ->setCategories("")->result;
my $categories = $soap -> setCategories ("GOTERM_BP_ALL,GOTERM_MF_ALL,GOTERM_CC_ALL,KEGG_PATHWAY");
if ($verbose_arg) {
  printf STDERR "II\tCategories are valid?:  $categories\n";  
}


########################################
##  Output table
########################################

my $chartReport = $soap -> getChartReport ($EASE_THRESHOLD, $GENE_COUNT);
my @chartRecords = $chartReport -> paramsout;

# shift(@chartRecords,($chartReport->result));
# print $chartReport->result."\n";

if ($verbose_arg) {
  printf STDERR "II\tTotal chart records:  %u\n\n", scalar (@chartRecords);
}

# my $retval = %{$chartReport->result};

# my @chartRecordKeys = keys %{$chartReport->result};
# my @chartRecordValues = values %{$chartReport->result};

print "Category\tTerm\tCount\t%\tPValue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";

my %chartRecord = %{$chartReport->result};
my $categoryName = $chartRecord{"categoryName"};
my $termName = $chartRecord{"termName"};
my $listHits = $chartRecord{"listHits"};
my $percent = $chartRecord{"percent"};
my $ease = $chartRecord{"ease"};
my $Genes = $chartRecord{"geneIds"};
my $listTotals = $chartRecord{"listTotals"};
my $popHits = $chartRecord{"popHits"};
my $popTotals = $chartRecord{"popTotals"};
my $foldEnrichment = $chartRecord{"foldEnrichment"};
my $bonferroni = $chartRecord{"bonferroni"};
my $benjamini = $chartRecord{"benjamini"};
my $FDR = $chartRecord{"afdr"};

##  Clean up the genes
$Genes = MapToNames ($Genes);
  
printf STDOUT "%s", $categoryName;
printf STDOUT "\t%s", $termName;
printf STDOUT "\t%u", $listHits;
printf STDOUT "\t%f", $percent;
printf STDOUT "\t%e", $ease;
printf STDOUT "\t%s", $Genes;
printf STDOUT "\t%u", $listTotals;
printf STDOUT "\t%u", $popHits;
printf STDOUT "\t%u", $popTotals;
printf STDOUT "\t%f", $foldEnrichment;
printf STDOUT "\t%e", $bonferroni;
printf STDOUT "\t%e", $benjamini;
printf STDOUT "\t%e", $FDR;
printf STDOUT "\n";

for (my $k = 0; $k < scalar (@chartRecords); $k++) {
  %chartRecord = %{$chartRecords[$k]};
  $categoryName = $chartRecord{"categoryName"};
  $termName = $chartRecord{"termName"};
  $listHits = $chartRecord{"listHits"};
  $percent = $chartRecord{"percent"};
  $ease = $chartRecord{"ease"};
  $Genes = $chartRecord{"geneIds"};
  $listTotals = $chartRecord{"listTotals"};
  $popHits = $chartRecord{"popHits"};
  $popTotals = $chartRecord{"popTotals"};
  $foldEnrichment = $chartRecord{"foldEnrichment"};
  $bonferroni = $chartRecord{"bonferroni"};
  $benjamini = $chartRecord{"benjamini"};
  $FDR = $chartRecord{"afdr"};
  
  ##  Clean up the genes
  $Genes = MapToNames ($Genes);
  
  printf STDOUT "%s", $categoryName;
  printf STDOUT "\t%s", $termName;
  printf STDOUT "\t%u", $listHits;
  printf STDOUT "\t%f", $percent;
  printf STDOUT "\t%e", $ease;
  printf STDOUT "\t%s", $Genes;
  printf STDOUT "\t%u", $listTotals;
  printf STDOUT "\t%u", $popHits;
  printf STDOUT "\t%u", $popTotals;
  printf STDOUT "\t%f", $foldEnrichment;
  printf STDOUT "\t%e", $bonferroni;
  printf STDOUT "\t%e", $benjamini;
  printf STDOUT "\t%e", $FDR;
  printf STDOUT "\n";
}        
  
