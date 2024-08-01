#!/home/charmm-gui/local/miniconda3/bin/python
# Identify non-aromatic rings.
# Copyright © 2022 Dan T. Major

import sys
import os
from rdkit import Chem
from openbabel import pybel

def write_header():
    print("\nRing finder program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

def write_footer():
    print("\nEnd ring finder program, returning to EnzyDock main\n")

def print_part1(file):
    file.write(
'''* Stream in SAMD boolean. Set to false if no flexible ring found.
*
\n''')

def print_part2(file):
    file.write('''\nreturn\n''')

def handle_input():
    # Given a pdb file, translate into Smiles string
    # and check whether non-aromatic rings are present
    pdb_file = sys.argv[1]
    pdb_file = pdb_file.lower()   # CHARMM sends in upper case parameters, change to lower
    print("Parsing ligand file: ",pdb_file)

    lig_id = sys.argv[2]
    #samdfile = "samd" + str(lig_id) + ".str"
    ## Write the samd_{lig_id}.str file for EnzyDock
    #file = open(samdfile,"w")
    #print_part1(file)
    #print("Writing file: ", samdfile)

    #try:
    bck_filename = correct_pdb(pdb_file)
    # Catch problems with molecule (RDKit returns None if there's a problem)
    mol_tmp = Chem.MolFromPDBFile(pdb_file)
    if mol_tmp != None:
       mol_pdb = mol_tmp
       num_atoms = mol_pdb.GetNumAtoms(onlyExplicit=True)
       smiles = Chem.MolToSmiles(mol_pdb,True,False)
       print("Ligand Smiles: ",smiles)
    else:
       print("Problems with molecule %s ..." % pdb_file)
    #mol_pdb = Chem.MolFromPDBFile(pdb_file)
    #num_atoms = mol_pdb.GetNumAtoms(onlyExplicit=True)
    #except:
       print("Error in RDKit read of ligand pdb file %s ..." % pdb_file)
       smiles = None 
       mol_pdb = None
       #Restore original file
       os.system('cp -f ' + bck_filename + ' ' + pdb_file)
       #print("Stopping EnzyDock ...\n")
       #file.write("\nstop\n")
       #print("Error in ring_aromatic.py", file=sys.stderr)
       #raise SystemExit(1)

    return smiles, mol_pdb, lig_id

def isRingAromatic(mol, bondRing):
    for id in bondRing:
       if not mol.GetBondWithIdx(id).GetIsAromatic():
          return False
    return True

def ring_aromatic(m):
    Chem.SanitizeMol(m)
    ri = m.GetRingInfo()
    # You can interrogate the RingInfo object to tell you the atoms that make up each ring:
    print(ri.AtomRings())
    # or the bonds that make up each ring:
    print(ri.BondRings())
    # To detect aromatic rings, loop over the bonds in each ring and
    # flag the ring as aromatic if all bonds are aromatic:
    flex_ring = False
    for i in range(len(ri.BondRings())):
        if (isRingAromatic(m, ri.BondRings()[i])):
            print("Ring %d is aromatic ..." % (i+1))
        else:
            print("Ring %d is not aromatic ..." % (i+1))
            flex_ring = True
    return flex_ring

def correct_pdb(filename):
    for mymol in pybel.readfile("pdb", filename):
        # First backup original pdb file, then correct pdb file
        exist = True
        i = 0
        while (exist):
           bck_filename = filename + '_' + str(i)
           if (not os.path.exists(bck_filename)):
              os.system('cp -f ' + filename + ' ' + bck_filename)
              exist = False
           else:
              i += 1
        output = pybel.Outputfile("pdb", filename, overwrite=True)
        print(output.write(mymol))
        output.close()
    return bck_filename

def write_samd_file(flex_ring, lig_id):
    samdfile = "../results/samd_" + str(lig_id) + ".str"
    # Write the samd_{lig_id}.str file for EnzyDock
    file = open(samdfile,"w")
    print("Writing file: ", samdfile)
    print_part1(file)
    line = "set samd = " + str(str(flex_ring).lower()) + "\n"
    file.write(line)
    print_part2(file)

def main():
    write_header()
    smiles, mol_pdb, lig_id = handle_input()
    # If we're unable to determine whether there's an unsaturated ring
    # (e.g., RDKit or pybel bug), just do SAMD
    if smiles != None:
       flex_ring = ring_aromatic(mol_pdb)
    else:
       flex_ring = True
    write_samd_file(flex_ring, lig_id)
    write_footer()

if __name__ == "__main__":
    main()

