rule Download_Combination:
  input:
    Expand_Combination
  output:
    output_fn1="Data/{combination}.download"
  run:
    from pathlib import Path

    ###  Read in the mapping of SLURM IDs to abbreviated sample names
    map_filename = {}
    for sample_name in config['samples']:
      slurm_id = config['samples'][sample_name]['fpkm']
      map_filename[slurm_id] = sample_name
    
    ###  Rename the SLURM IDs to the abbreviated sample names only if the source exists and the target doesn't already exist
    for key in map_filename.keys ():
      if os.path.exists ("Data/fpkm/{src}".format (src=key)):
        if not os.path.exists ("Data/fpkm/{target}".format (target=map_filename[key])):
          os.rename ("Data/fpkm/{src}".format (src=key), "Data/fpkm/{target}".format (target=map_filename[key]))
          
    Path (output.output_fn1).touch ()
    

rule Download_Sample:
  input:
    input_fn1="Data/status/{sample}.ballgown",
    input_fn2="Data/status/{sample}.mapping",
    input_fn3="Data/status/{sample}.summary"
  output:
    output_fn1="Data/status/{sample}.all"
  run:
    from pathlib import Path
    
    Path (output.output_fn1).touch ()


rule Download_Ballgown:
  input:
    ReturnEmpty
  output:
    output_fn1="Data/status/{sample}.ballgown"
  params:
    sample="{sample}",
    fpkm_slurm_id = lambda wildcards: config["samples"][wildcards.sample]["fpkm"]
  run:
    import os
    from pathlib import Path

    if not os.path.exists ("Data/"):
      os.makedirs ("Data/")

    if not os.path.exists ("Data/fpkm/"):
      os.makedirs ("Data/fpkm/")

    shell ("{cmd} {login}@{server}:{path}/{id}/limit_ballgown/* Data/fpkm/{s}/".format (cmd=RSYNC_SO, login=DATA_LOGIN, server=DATA_SERVER, path=DATA_PATH, id=params.fpkm_slurm_id, s=params.sample))
    
    Path (output.output_fn1).touch ()
    

rule Download_Mapping:
  input:
    ReturnEmpty
  output:
    output_fn1="Data/status/{sample}.mapping"
  params:
    sample="{sample}",
    fpkm_slurm_id = lambda wildcards: config["samples"][wildcards.sample]["fpkm"],
    read_prefix = lambda wildcards: config["samples"][wildcards.sample]["readprefix"]
  run:
    import os
    from pathlib import Path

    if not os.path.exists ("Data/"):
      os.makedirs ("Data/")

    if not os.path.exists ("Data/mapping/"):
      os.makedirs ("Data/mapping/")
      
    shell ("{cmd} {login}@{server}:{path}/{id}/{srr}.flagstat Data/mapping/{s}/".format (cmd=RSYNC_SO, login=DATA_LOGIN, server=DATA_SERVER, path=DATA_PATH, id=params.fpkm_slurm_id, srr=params.read_prefix, s=params.sample))
    
    Path (output.output_fn1).touch ()


rule Download_Summary:
  input:
    ReturnEmpty
  output:
    output_fn1="Data/status/{sample}.summary"
  params:
    sample="{sample}",
    summary_slurm_id = lambda wildcards: config["samples"][wildcards.sample]["summary"],
    read_prefix = lambda wildcards: config["samples"][wildcards.sample]["readprefix"]
  run:
    import os
    from pathlib import Path

    if not os.path.exists ("Data/"):
      os.makedirs ("Data/")

    if not os.path.exists ("Data/summary/"):
      os.makedirs ("Data/summary/")
      
    shell ("{cmd} {login}@{server}:{path}/{id}/{srr}*.txt Data/summary/{s}/".format (cmd=RSYNC_SO, login=DATA_LOGIN, server=DATA_SERVER, path=DATA_PATH, id=params.summary_slurm_id, srr=params.read_prefix, s=params.sample))
    
    Path (output.output_fn1).touch ()


