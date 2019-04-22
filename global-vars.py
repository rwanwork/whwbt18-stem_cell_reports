#####################################################################
##  Global variables
##
##  Raymond Wan (raymondwan@ust.hk)
##  Organizations:
##    - Division of Life Science, 
##      Hong Kong University of Science and Technology
##      Hong Kong
##
##  Copyright (C) 2018, Raymond Wan, All rights reserved.
#####################################################################

##  Login information for the server where the data is stored
DATA_LOGIN = ""
DATA_SERVER = ""
DATA_PATH = ""

##  Login information for DAVID
EMAIL = ""


#####################################################################
##  Path to programs
#####################################################################

RSYNC_SO = "rsync --msgs2stderr --cvs-exclude --itemize-changes --progress -azv --size-only"


#####################################################################
##  Data samples to process
#####################################################################

##  Samples table
SAMPLES_DIRECTORIES = ["meta", "mapping", "summary"]

