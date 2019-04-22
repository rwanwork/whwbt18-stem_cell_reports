rule Symlink:
  input:
    input_fn1="Data/{combination}.download"
  output:
    output_fn1="Data/{combination}.symlink",
    output_fn2="Data/combos/{combination}.csv"
  params:
    combination="{combination}"
  run:
    import os
    from pathlib import Path

    target = "Data/" + params.combination + "/"
    if not os.path.exists (target):
      os.makedirs (target)

    for s in config['combinations'][params.combination]:
      src = "../fpkm/" + s
      target = "Data/" + params.combination + "/" + s
      os.symlink (src, target)

    text_dir = "Text/" + params.combination + "/ballgown/"
    if not os.path.exists (text_dir):
      os.makedirs (text_dir)
    
    graphs_dir = "Graphs/" + params.combination + "/ballgown/"
    if not os.path.exists (graphs_dir):
      os.makedirs (graphs_dir)

    ##  Create the CSV file that is required by ballgown
    if not os.path.exists ("Data/combos/"):
      os.makedirs ("Data/combos/")
      
    fp = open (output.output_fn2, 'w')
    fp.write ("\"ids\",\"condition\"" + "\n")
    for s in config['combinations'][params.combination]:
      out_str = "\"" + s + "\",\"" + config['samples'][s]['condition'] + "\"\n"
      fp.write (out_str)
    fp.close ()
    
    Path (output.output_fn1).touch ()
    

rule Ballgown:
  input:
    input_fn1="Data/{combination}.symlink",
    input_fn2="Data/combos/{combination}.csv"
  output:
    output_fn1="Text/{combination}/ballgown/genes-all.tsv",
    output_fn2="Text/{combination}/ballgown/genes-mean.tsv",
    output_fn3="Text/{combination}/ballgown/genes-fpkm.tsv"
  params:
    combination="{combination}"
  shell:
    """
    source activate pompe
    
    R/ballgown.R --sample {params.combination} --results Text/{params.combination}/ballgown/ --graphs Graphs/{params.combination}/ballgown/
    
    cd Text/{params.combination}/ballgown/
    cd ../../../

    ##  Remove the CSV file that is required by ballgown
    rm -f {input.input_fn2}
    """


rule Hierarchical_Clustering:
  input: 
    input_fn="Text/{combination}/ballgown/genes-fpkm.tsv"
  output:
    output_fn="Graphs/{combination}/HC/hc-pearson.eps"
  params:
    combination="{combination}"
  shell:
    """
    source activate pompe

    R/hierarchical-clustering.R --input {input.input_fn} --graphs Graphs/{params.combination}/HC/
    """

  
rule PCA:
  input: 
    input_fn="Text/{combination}/ballgown/genes-fpkm.tsv"
  output:
    output_fn="Graphs/{combination}/PCA/pca.eps"
  params:
    combination="{combination}"
  shell:
    """
    source activate pompe

    R/pca.R --input {input.input_fn} --graphs Graphs/{params.combination}/PCA/
    """

  
rule tSNE:
  input: 
    input_fn="Text/{combination}/ballgown/genes-fpkm.tsv"
  output:
    output_fn="Graphs/{combination}/tSNE/tsne.eps"
  params:
    combination="{combination}"
  shell:
    """
    source activate pompe

    R/tsne.R --input {input.input_fn} --graphs Graphs/{params.combination}/tSNE/
    """
  

