#####################################################################
##  Snakefile
##
##  Raymond Wan (raymondwan@ust.hk)
##  Organizations:
##    - Division of Life Science, 
##      Hong Kong University of Science and Technology
##      Hong Kong
##
##  Copyright (C) 2017, Raymond Wan, All rights reserved.
#####################################################################

##  Configuration file
configfile: "config.yaml"

##  Define global constraints on wildcards
wildcard_constraints:
  slurm_id = "\d+",
  combination = "combo\d"

##  Include additional functions and rules
include: "global-vars.py"
include: "general.py"  ##  General functions that may be useful in multiple locations
include: "expand.py"
include: "samples-primitives.py"
include: "samples.py"  ##  Functions and rules for creating a table of samples
include: "download.py"
include: "ballgown.py"
include: "pairwise.py"
  
  
##  Beginning of rules
rule all:
  input: 
    "Complete/combo3.samples",
    "Complete/combo7.samples",
    "Complete/combo3/MPC-ASC.done",
    "Complete/combo7/MPC-ASC.done",
    "Complete/combo3.done",
    "Complete/combo7.done"


rule Complete_Pairwise:
  input:
    input_fn1="Graphs/{combination}/{sample1}-{sample2}/david/kegg.eps",
    input_fn2="Graphs/{combination}/{sample1}-{sample2}/david/go.eps",
    input_fn3="Graphs/{combination}/{sample1}-{sample2}/maplot/maplot.eps",
    input_fn4="Graphs/{combination}/{sample1}-{sample2}/volcano/volcano.eps"    
  output:
    output_fn1="Complete/{combination}/{sample1}-{sample2}.done"
  shell:
    """
    touch {output.output_fn1}
    """
    

rule Complete_Ballgown:
  input:
    input_fn1="Graphs/{combination}/HC/hc-pearson.eps",
    input_fn2="Graphs/{combination}/PCA/pca.eps"
    #input_fn3="Graphs/{combination}/tSNE/tsne.eps"
  output:
    output_fn1="Complete/{combination}.done"
  shell:
    """
    touch {output.output_fn1}
    """
    
  
rule Complete_SamplesList:
  input: 
    input_fn1="Text/{combination}/samples/data-table.csv"
  output:
    output_fn="Complete/{combination}.samples"
  shell:
    """
    touch {output.output_fn}
    """


