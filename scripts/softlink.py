#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
# Create softlinks to allow long file names with CHARMM
# Copyright © 2022 Dan T. Major

import sys
import os
from os import path

print("\nSoft-link program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

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

print("\nEnd soft-link program, returning to EnzyDock main\n")

