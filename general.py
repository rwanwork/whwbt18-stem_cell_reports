#####################################################################
##  General Python functions
##
##  Raymond Wan (raymondwan@ust.hk)
##  Organizations:
##    - Division of Life Science, 
##      Hong Kong University of Science and Technology
##      Hong Kong
##
##  Copyright (C) 2018, Raymond Wan, All rights reserved.
#####################################################################

##  Convert a integral number to "millions"
def Millions (x):
  if x != "N/A":
    x = str (round (int (x) / 1000 / 1000, 1))
  
  return (x)  


##  Convert a integral number to "billions"
def Billions (x):
  if x != "N/A":
    x = str (round (int (x) / 1000 / 1000 / 1000, 1))
  
  return (x)  

