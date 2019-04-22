def Expand_Combination (wc):
  all_paths = []
  for sample in config['combinations'][wc.combination]:
    #slurm_id = config['samples'][sample]
    #path = "Data/status/" + str (slurm_id)
    path = "Data/status/" + sample + ".all"
    all_paths.append (path)
  
  print ("Expand_Combinations:\t", all_paths, file=sys.stderr)
  
  return all_paths

    
def ReturnEmpty (wc):
  return []  


