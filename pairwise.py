##  Performing averaging over replicates.  Irrelevant columns are removed and a t-test is performed.
rule Average_Replicates:
  input: 
    input_fn="Text/{combination}/ballgown/genes-fpkm.tsv"
  output:
    output_fn1="Text/{combination}/{sample1}-{sample2}/significant/genes-sig-all.tsv",
    output_fn2="Text/{combination}/{sample1}-{sample2}/significant/genes-sig-filtered.tsv"
  params:
    sample1="{sample1}",
    sample2="{sample2}"
  shell:
    """
    source activate pompe
    
    R/avg-replicates.R --input {input.input_fn} --allgenes {output.output_fn1} --siggenes {output.output_fn2} --sample1 {params.sample1} --sample2 {params.sample2} 
    """


rule MAPlot:
  input: 
    input_fn1="Text/{combination}/{sample1}-{sample2}/significant/genes-sig-all.tsv"
  output:
    output_fn1="Graphs/{combination}/{sample1}-{sample2}/maplot/maplot.eps"
  params:
    combination="{combination}",
    sample1="{sample1}",
    sample2="{sample2}"
  shell:
    """
    source activate pompe
    
    R/maplot.R --input {input.input_fn1} --output {output.output_fn1} --sample1 {params.sample1} --sample2 {params.sample2}
    """


rule Volcano:
  input: 
    input_fn1="Text/{combination}/{sample1}-{sample2}/significant/genes-sig-all.tsv"
  output:
    output_fn1="Graphs/{combination}/{sample1}-{sample2}/volcano/volcano.eps"
  params:
    combination="{combination}",
    sample1="{sample1}",
    sample2="{sample2}"
  shell:
    """
    source activate pompe
    
    R/volcano-plot.R --input {input.input_fn1} --output {output.output_fn1} --sample1 {params.sample1} --sample2 {params.sample2}
    """

    
rule Calculate_DAVID:
  input:
    input_fn="Text/{combination}/{sample1}-{sample2}/significant/genes-sig-filtered.tsv"
  output:
    output_fn1="Text/{combination}/{sample1}-{sample2}/david/david-up.tsv",
    output_fn2="Text/{combination}/{sample1}-{sample2}/david/david-down.tsv",
    output_fn3="Text/{combination}/{sample1}-{sample2}/david/david-combined.tsv"
  shell:
    """
    source activate pompe
    
    perl Perl/enrichment-david.pl --email {EMAIL} --regulate up --verbose <{input.input_fn} >{output.output_fn1} 2>{output.output_fn1}.stderr
    sleep 5m

    perl Perl/enrichment-david.pl --email {EMAIL} --regulate down --verbose <{input.input_fn} >{output.output_fn2} 2>{output.output_fn2}.stderr
    sleep 5m
    
    perl Perl/enrichment-david.pl --email {EMAIL} --regulate combined --verbose <{input.input_fn} >{output.output_fn3} 2>{output.output_fn3}.stderr
    sleep 5m
    """


rule Plot_DAVID:
  input: 
    input_fn1="Text/{combination}/{sample1}-{sample2}/david/david-up.tsv",
    input_fn2="Text/{combination}/{sample1}-{sample2}/david/david-down.tsv",
    input_fn3="Text/{combination}/{sample1}-{sample2}/david/david-combined.tsv"
  output:
    output_fn1="Graphs/{combination}/{sample1}-{sample2}/david/kegg.eps",
    output_fn2="Graphs/{combination}/{sample1}-{sample2}/david/go.eps",
    output_fn3="Graphs/{combination}/{sample1}-{sample2}/david/kegg-combined.eps",
    output_fn4="Graphs/{combination}/{sample1}-{sample2}/david/go-combined.eps"
  shell:
    """
    source activate pompe
    
    R/funcchart.R --upregulate {input.input_fn1} --downregulate {input.input_fn2} --kegg {output.output_fn1} --go {output.output_fn2}
    
    ##  If we're in "combined" mode, then just provide the same arguments to
    ##    --upregulate and --downregulate
    R/funcchart.R --upregulate {input.input_fn3} --downregulate {input.input_fn3} --kegg {output.output_fn3} --go {output.output_fn4}
    """


