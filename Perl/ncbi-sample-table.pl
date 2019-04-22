#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;

my %names;
my %name_to_readprefix;
my %readprefix_to_name;
my %readprefix_to_submission;
my %readprefix_to_experiment;
my %readprefix_to_sample;
my %readprefix_to_study;
my %readprefix_to_biosample;
my %readprefix_to_bioproject;
my %readprefix_to_replacedby;

my %experiments_lookup;
my %experiment_to_gsm;

my %submissions_lookup;
my %submission_to_gse;

##  Search for the list of samples
while (<STDIN>) {
  my $line = $_;
  if ($line =~ /^samples:/) {
    last;
  }
}

while (<STDIN>) {
  my $sample = $_;
  my $fpkm = <STDIN>;
  my $summary = <STDIN>;
  my $readprefix = <STDIN>;
  my $condition = <STDIN>;
  if (!defined ($fpkm)) {
    last;
  }

  chomp ($sample);
  chomp ($fpkm);
  chomp ($summary);
  chomp ($readprefix);
  chomp ($condition);
  
  if ($sample =~ /\s+(.+):$/) {
    $sample = $1;
  }
  
  if ($readprefix =~ /\s+readprefix:  (.+)$/) {
    $readprefix = $1;
  }

  if (defined $names{$sample}) {
    printf STDERR "EE\tRepeated sample found:  %s\n", $sample;
    exit (1);
  }
  
  $names{$sample} = $sample;
  $name_to_readprefix{$sample} = $readprefix;
  $readprefix_to_name{$readprefix} = $sample;
}


##  First pass -- get all the information that we can get using the readprefix and type RUN
my $fn = "./SRA_Accessions.tab";
open (my $fp, "<", $fn) or die "EE\tCould not find $fn.\n";
while (<$fp>) {
  my $line = $_;
  chomp ($line);
  
  my ($accession, $submission, $status, $updated, $published, $received, $type, $center, $visibility, $alias, $experiment, $sample, $study, $loaded, $spots, $bases, $md5sum, $biosample, $bioproject, $replacedby) = split /\t/, $line;
  if ($type ne "RUN") {
    next;
  }

  if (defined ($readprefix_to_name{$accession})) {
    my $key = $accession;
    
    $submissions_lookup{$submission} = $submission;
    $experiments_lookup{$experiment} = $experiment;

    $readprefix_to_submission{$key} = $submission;
    $readprefix_to_experiment{$key} = $experiment;
    $readprefix_to_sample{$key} = $sample;
    $readprefix_to_study{$key} = $study;
    $readprefix_to_biosample{$key} = $biosample;
    $readprefix_to_bioproject{$key} = $bioproject;
    $readprefix_to_replacedby{$key} = $replacedby;
  }
}
close ($fp);


##  Second pass -- get all the information that we can using the type EXPERIMENT
open ($fp, "<", $fn) or die "EE\tCould not find $fn.\n";
while (<$fp>) {
  my $line = $_;
  chomp ($line);
  
  my ($accession, $submission, $status, $updated, $published, $received, $type, $center, $visibility, $alias, $experiment, $sample, $study, $loaded, $spots, $bases, $md5sum, $biosample, $bioproject, $replacedby) = split /\t/, $line;
  if ($type ne "EXPERIMENT") {
    next;
  }

  if (defined ($experiments_lookup{$accession})) {
    my $key = $accession;
    
    $experiment_to_gsm{$key} = $alias;
  }
}
close ($fp);


##  Second pass -- get all the information that we can using the type SUBMISSION
open ($fp, "<", $fn) or die "EE\tCould not find $fn.\n";
while (<$fp>) {
  my $line = $_;
  chomp ($line);
  
  my ($accession, $submission, $status, $updated, $published, $received, $type, $center, $visibility, $alias, $experiment, $sample, $study, $loaded, $spots, $bases, $md5sum, $biosample, $bioproject, $replacedby) = split /\t/, $line;
  if ($type ne "SUBMISSION") {
    next;
  }

  if (defined ($submissions_lookup{$accession})) {
    my $key = $accession;
    
    if ($alias =~ /^GEO: (.+)$/) {
      $alias = $1;
    }
    
    $submission_to_gse{$key} = $alias;
  }
}
close ($fp);

foreach my $tmp_key (sort (keys %names)) {
  printf "%s", $names{$tmp_key};
  if ($names{$tmp_key} =~ /^MPC_/) {
    printf "\n";
    next;
  }

  ##  Change key to the SRR number
  my $key = $name_to_readprefix{$tmp_key};
  printf "\t%s", $key;
  printf "\t%s", $readprefix_to_experiment{$key};
  printf "\t%s", $experiment_to_gsm{$readprefix_to_experiment{$key}};
  printf "\t%s", $readprefix_to_submission{$key};
  printf "\t%s", $submission_to_gse{$readprefix_to_submission{$key}};
  printf "\t%s", $readprefix_to_sample{$key};
  printf "\t%s", $readprefix_to_study{$key};
  printf "\t%s", $readprefix_to_biosample{$key};
  printf "\t%s", $readprefix_to_bioproject{$key};
  printf "\t%s", $readprefix_to_replacedby{$key};
  
  printf "\n";
}


