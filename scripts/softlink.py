#!/home/qnt/majort/anaconda3/bin/python3.6

import sys
import os
from os import path

# Handle inpit arguments
input_arg = sys.argv
arg1 = input_arg[1]#.lower()
arg2 = input_arg[2].lower().rstrip('/')   # CHARMM sends in upper case parameters, change to lower
arg3 = int(input_arg[3])

# Create softlink at beginning of job and remove after
if (arg3 == 1):
   # First remove link
   if path.islink(arg2):
      command = "unlink " + arg2
      print(command)
      os.system(command)
   # Command to  create link
   command = "ln -s " + arg1 + "/" + arg2 + " " + arg2
elif (arg3 == 0):
   command = "unlink " + arg2

print(command)
os.system(command)

