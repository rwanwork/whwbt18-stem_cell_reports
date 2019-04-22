rule Calculate_Samples_Metadata:
  input:
    input_flag="Data/status/{sample}.summary"
  output: 
    output_fn="Text/{combination}/samples/meta/{sample}.txt"
  run:
    import sys

    fn = output.output_fn
    out_fp = open (fn, 'w')
    out_fp.write (wildcards.sample + "\n")
    out_fp.close ()


rule Calculate_Samples_Mapping:
  input: 
    input_flag="Data/status/{sample}.summary"
  output: 
    output_fn="Text/{combination}/samples/mapping/{sample}.txt"
  params:
    sample="{sample}",
    read_prefix = lambda wildcards: config["samples"][wildcards.sample]["readprefix"]
  run:
    import re
    
    input_fn = "Data/mapping/" + params.sample + "/" + params.read_prefix + ".flagstat"
    in_fp = open (input_fn, 'r')
    empty_string = ''

    line = in_fp.readline ()  ##  Read and discard the header line
    line = in_fp.readline ()

    mapping_rate = ""
    #  "properly paired" doesn't work since there are single-ended reads    
#    pattern = re.compile ('properly paired \((?P<rate>.+)%')
    pattern = re.compile ('mapped \((?P<rate>.+)%')
    while line != empty_string:
      line = line.rstrip ()
  
      if pattern.search (line):
        result = pattern.search (line)
        mapping_rate = result.group ('rate')
        #print (mapping_rate, "\n")
        break
      line = in_fp.readline ()

    in_fp.close ()

    ##  Write out the value    
    fn = output.output_fn
    out_fp = open (fn, 'w')
    out_fp.write (mapping_rate + "\n")
    out_fp.close ()
  

rule Calculate_Samples_Summary:
  input: 
    input_flag="Data/status/{sample}.summary"
  output: 
    output_fn="Text/{combination}/samples/summary/{sample}.txt"
  params:
    sample="{sample}",
    read_prefix = lambda wildcards: config["samples"][wildcards.sample]["readprefix"]
  run:
    import re
    import os.path

    input_fn1 = "Data/summary/" + params.sample + "/" + params.read_prefix + "_1.txt"
    input_fn2 = "Data/summary/" + params.sample + "/" + params.read_prefix + "_2.txt"

    input_fns = []
    input_fns.append (input_fn1)
    input_fns.append (input_fn2)

    if os.path.exists (input_fn2):
      out_str = "Paired-end"
    else:
      out_str = "Single-end"

    for fn in input_fns:
      empty_string = ''
      
      numreads = "N/A"
      minlen = "N/A"
      maxlen = "N/A"
      avglen = "N/A"
      numbases = "N/A"
    
      numreads_pattern = re.compile ('Number of reads:\s+(?P<numreads>.+)')
      minlen_pattern = re.compile ('Minimum read length:\s+(?P<minlen>.+)')
      maxlen_pattern = re.compile ('Maximum read length:\s+(?P<maxlen>.+)')
      avglen_pattern = re.compile ('Average read length:\s+(?P<avglen>.+)')
      numbases_pattern = re.compile ('Total number of bases:\s+(?P<numbases>.+)')
    
      if os.path.exists (fn):
        in_fp = open (fn, 'r')

        line = in_fp.readline ()  ##  Read and discard the header line
        line = in_fp.readline ()

        while line != empty_string:
          line = line.rstrip ()
  
          if numreads_pattern.search (line):
            result = numreads_pattern.search (line)
            numreads = result.group ('numreads')
          if minlen_pattern.search (line):
            result = minlen_pattern.search (line)
            minlen = result.group ('minlen')
          if maxlen_pattern.search (line):
            result = maxlen_pattern.search (line)
            maxlen = result.group ('maxlen')
          if avglen_pattern.search (line):
            result = avglen_pattern.search (line)
            avglen = result.group ('avglen')
          if numbases_pattern.search (line):
            result = numbases_pattern.search (line)
            numbases = result.group ('numbases')
          line = in_fp.readline ()
          
        in_fp.close ()

      out_str = out_str + "\t" + Millions (numreads) + "\t" + minlen + "\t" + maxlen + "\t" + avglen + "\t" + Billions (numbases);

    ##  Write out the value    
    fn = output.output_fn
    out_fp = open (fn, 'w')
    out_fp.write (out_str + "\n")
    out_fp.close ()


rule Generate_Samples_Record:
  input:
    Create_Samples_Record
  output:
    output_fn="Text/{combination}/samples/records/{sample}.txt"
  run:
    directories = SAMPLES_DIRECTORIES
    
    ##  Where to write the output file
    out_fn = output.output_fn
    out_fp = open (out_fn, 'w')
    
    base_string = "Text/" + wildcards.combination + "/samples/{0}/{1}.txt"
    all_str = ""
    for curr_data in directories:
      curr_fn = base_string.format (curr_data, wildcards.sample)
      #print ("[", curr_fn, "]\n", file=sys.stderr)

      curr_fp = open (curr_fn, 'r')
      curr_str = curr_fp.read ();
      curr_str = curr_str.rstrip ();  ##  Remove trailing newline
      all_str = all_str + "\t" + curr_str
      curr_fp.close ()

    all_str = all_str.lstrip ()  ##  Remove the newline at the front
    out_fp.write (all_str + "\n")
    
    ##  Close the output file when done
    out_fp.close ()
    

##  Create the entire data table
rule Generate_Samples_CSV:
  input: 
    Expand_Samples_CSV
  output:
    output_fn="Text/{combination}/samples/data-table.csv"
  shell:
    """
    perl Perl/create-samples-header.pl >{output}    
    cat {input} >>{output.output_fn}
    """
    
