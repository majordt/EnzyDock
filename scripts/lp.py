#!/home/charmm-gui/local/miniconda3/bin/python

# identify lone pair in cgenff topology file and omit it.
# Transfer charge to bonding atom 
# call it from enzydock the first time to ligand is read.

import subprocess
import sys

def handle_input():
    # Handle input arguments and return values
    input_arg = sys.argv
    lig = str(input_arg[1]) # @currligand
    return lig

def lp_name(lig):
  ligf = "../local_top/ligand_"+str(lig)+".str"
  with open(ligf) as f:
    ls = f.readlines()
  lps = []
  halogens = []
  for l in ls:
    if l.upper().startswith("LONEPAIR"):
      lps.append(l.split()[2]) #atom name of lone pair
      halogens.append(l.split()[3]) #atom name of halogen      
  return lps, halogens

def get_charges(lps, halogens, lig):
  ligf = "../local_top/ligand_"+str(lig)+".str"
  with open(ligf) as f:
    ls = f.readlines()
  charges = {}
  for l in ls:
    if l.upper().startswith("ATOM"):
      atom = l.split()[1]
      charge = l.split()[3]
      if atom in (lps+halogens):
        charges[atom] = charge
  return charges

def merge_charges(lps, halogens, charges):
  new_charges = {}
  for atom,lp in zip(halogens,lps):
    new_charges[atom] = str(round(float(charges[atom]) + float(charges[lp]),3))
  return new_charges
  
def re_write_str(lps, halogens, new_charges,lig):
  ligf = "../local_top/ligand_"+str(lig)+".str"
  with open(ligf) as f:
    ls = f.readlines()
  new_lines = ""
  for l in ls:
    new_line = l
    if l.startswith("ATOM"):
      if l.split()[1] in lps: #comment out lone pairs
        new_line = "!"+l
      elif l.split()[1] in halogens: #adjust charge of halogens
        new_line = l.split()
        new_line[3] = new_charges[l.split()[1]] #chrage,atom name
        new_line = " ".join(new_line)+"\n"
    elif l.startswith("LONEPAIR"):
      new_line = "!"+l
    new_lines += new_line

  ligf0 = "../local_top/ligand_"+str(lig)+"_lp.str"
  subprocess.call(["/bin/mv",ligf,ligf0])

  with open(ligf, 'w') as f:
    f.write(new_lines)

def fix_pdb(lig, lps):
  # remove lonepai line from pdb file - assume it's the last of all atoms so
  #   connect section isn't messed up
  lps_l = [i.lower() for i in lps]
  ligf = "../pdb/ligand_"+str(lig)+".pdb"
  with open(ligf) as f:
    ls = f.readlines()
  new_lines = ""
  for l in ls:
    if l.startswith("ATOM"):
      if l[12:16].strip().lower() in lps_l: #delete lone pairs
        continue
      elif 'lp' in l[12:16].strip().lower():
        continue
    new_lines += l

  ligf0 = "../pdb/ligand_"+str(lig)+"_lp.pdb"
  subprocess.call(["/bin/mv",ligf,ligf0])
  with open(ligf, 'w') as f:
    f.write(new_lines)

def main():
  lig = handle_input()
  lps, halogens = lp_name(lig)
  if len(lps) > 0:
    fix_pdb(lig,lps)
    charges = get_charges(lps, halogens, lig)
    new_charges = merge_charges(lps, halogens, charges)
    re_write_str(lps, halogens, new_charges,lig)

if __name__ == "__main__":
  main()
