def Create_Samples_Record (wc):
  files = []
  
  directories = SAMPLES_DIRECTORIES
  
  base_string = "Text/{0}/samples/{1}/{2}.txt"
  for curr_data in directories:
    files.append (base_string.format (wc.combination, curr_data, wc.sample))
          
  print ("Create_Samples_Record:\t", files, file=sys.stderr)
  
  return files


def Expand_Samples_CSV (wc):
  files = []

  base_string = "Text/" + wc.combination + "/samples/records/{0}.txt"
  for s in config['combinations'][wc.combination]:
    fn = base_string.format (s)
    files.append (fn)
  
  print ("Expand_Samples_Record:\t", files, file=sys.stderr)
  
  return files


