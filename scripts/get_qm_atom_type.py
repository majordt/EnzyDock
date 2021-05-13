#!/home/qnt/majort/anaconda3/bin/python3.6

import sys
import os

def get_qm_atom_type():
    input_arg = sys.argv
    arg1 = input_arg[1].lower()   # CHARMM sends in upper case parameters, change to lower
    arg2 = input_arg[2].upper()   # Make sure it's upper
    file = open(arg1, 'r')
    lines = file.readlines()
    status = 0
    atom_type = ""
    for line in lines:
       if (len(line.strip()) > 0):
          line0 = line.split()
          if (line0[0].upper() == 'ATOM'):
             if (line0[1].upper() == arg2):
                atom_type = line0[2]
                status = 1
    return(atom_type, status)

def write_str_file(str_file, atom_type):
    str_file.write(
'''* QM boundary atom type
*

set qmligl ''')
    print(atom_type, file = str_file)

def main():
    atom_type, status = get_qm_atom_type()
    if os.path.exists('../local_top/qm_atom_type.str'):
       os.remove('../local_top/qm_atom_type.str')
    if (status == 1):
       file = open('../local_top/qm_atom_type.str', 'w')
       write_str_file(file, atom_type)
       file.close()

if __name__ == "__main__":
    main()

