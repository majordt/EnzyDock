#!/home/charmm-gui/local/miniconda3/bin/python
# Copyright © 2022 Dan T. Major

import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem import rdqueries
from rdkit.Chem import rdMolTransforms

def print_part1(file):
    file.write(
'''* Monte Carlo add moves
* Moves added must match the ones deleted in mc_delete.str
* Format is same as 'move' definition for covalent docking, with appropriate atom name modifications
* , increasing label names: tr2, tr3, ...
* Moves with translations ('rtrn') and rotations ('rrot') can also be defined, see commented lines. 
* Also, these may be combined into global moves (see 'move link' below).
*

! Weights are chosen such that each degree of freedom has roughly equal 
! frequency of being chosen.
! Rigid body translations of ligand
!move add mvtp rtrn byresidue weight 1.00 dmax 0.10 label atrn sele segid ligand_@currligand end 

! Rigid body rotations of ligand
!move add mvtp rrot byresidue weight 1.00 dmax 25.0 label arot sele atom ligand_@currligand 1 C25 end 

! Torsional rotations in ligand. Lines below first must be modified manually.
! Note that in the torsional definition there is DIRECTIONALITY!! The MC module will move the atoms after 
! the second atom in the definition. The MC module doesn't care about which atoms are fixed.
! So, if we have a covalently docked ligand via a Cys-S-C-…-X-Y-Z-W-… link, all MC moves MUST be defined 
! "outwards" and "away" from the Cys residue. Otherwise, the MC module will move the Cys residue itself!!!
! move add mvtp tors weight 1.00 dmax 180.0 label tr2  fewer 0 sele atom ligand_@currligand 1 Y  show end -
!                                                              sele atom ligand_@currligand 1 Z  show end
set ntors 0

if covalent .eq. true then
! move add mvtp tors weight 2.00 dmax 180.0 label tr1  fewer 0 sele atom FLX@{covseg} 1 @protlink0 show end -
!                                                              sele atom FLX@{covseg} 1 @protlink1 show end
! move add mvtp tors weight 2.00 dmax 180.0 label tr2  fewer 0 sele atom FLX@{covseg} 1 @protlink1 show end -

 move add mvtp tors weight 2.00 dmax 180.0 label tr1  fewer 0 sele atom FLX@{covseg} 1 @protlink1 show end -
                                                              sele atom ligand_@currligand 1 @liglink1 show end
! incr ntors by 2
 incr ntors by 1
endif

! Additional torsional definitions below
!! move add mvtp tors weight 1.00 dmax 180.0 label tr3 fewer 0 sele ... show end sele ... show end\n
! move add mvtp tors weight 1.00 dmax 180.0 label tr2 fewer 0 sele ... show end sele ... show end\n''')

def print_part2(file):
    file.write(
'''! Link moves so they are done together.
!move link lab1 tr1 lab2 tr2 

return\n''')

def print_part3(file):
    file.write(
'''* Monte Carlo delete moves
* Moves deleted must match the ones added in mc_add.str
* All moves must be unlinked before any delete takes place
*
\n''')

def print_part4(file):
    file.write(
'''if covalent .eq. true then
 move dele label tr1
! move dele label tr2
! decr ntors by 2
 decr ntors by 1
endif

! Additional torsional deletions below\n''')

def print_part5(file):
    file.write(
'''\nreturn\n''')

def write_header():
    print("\nTorsional finder program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

def write_footer():
    print("\nEnd torsional finder program, returning to EnzyDock main\n")

def smiles_tors():
    # Given a pdb file, translate into Smiles string (with possible given starting atom)
    # and generate a list of rotatable bonds. Rotatable bond list is written to CHARMM
    # stream files (mc_add.str, mc_delete.str) that are used by the MC module used by EnzyDock.
    # A crucial assumption is that the rotatable bonds are written in an "outwards" manner
    # when covalent docking is done (otherwise CHARMM could rotate protein and not ligand!)

    pdb_file = sys.argv[1]
    pdb_file = pdb_file.lower()   # CHARMM sends in upper case parameters, change to lower
    print("Parsing ligand file: ",pdb_file)

    mcfile1 = "../stream/mc/mc_add.str"
    mcfile2 = "../stream/mc/mc_delete.str"
    # Write the mc_add.str file for EnzyDock
    file = open(mcfile1,"w")
    print_part1(file)
    print("Writing file: ", mcfile1)

    #try:
    # Catch problems with molecule (RDKit returns None if there's a problem)
    mol_tmp = Chem.MolFromPDBFile(pdb_file)
    if mol_tmp != None:
       mol_pdb = mol_tmp
       num_atoms = mol_pdb.GetNumAtoms(onlyExplicit=True)
    else:
       print("Problems with molecule %s ..." % pdb_file)
    #mol_pdb = Chem.MolFromPDBFile(pdb_file)
    #num_atoms = mol_pdb.GetNumAtoms(onlyExplicit=True)
    #except:
       print("Error in RDKit read of ligand pdb file %s ..." % pdb_file)
       print("Stopping EnzyDock ...\n")
       file.write("\nstop\n")
       with open("enzydock.log",'a') as f:
         f.write("Problems with molecule %s ...\n" % pdb_file)
         f.write("Error in RDKit read of ligand pdb file %s ...\n" % pdb_file)
         f.write("Stopping EnzyDock ...\n")
       raise SystemExit(1)

    print("Number of ligand atoms: ",num_atoms)
    idx = -1
    ncov = 1   # Number of covalent torsions (Check with RS)
    # Find atom that serves as seed (covalently bonded)
    if len(sys.argv) > 2:
        for i in range(num_atoms):
            atom = mol_pdb.GetAtomWithIdx(i)
#        print(atom.GetPDBResidueInfo().GetName().strip())
            if (atom.GetPDBResidueInfo().GetName().strip() == sys.argv[2]):
                idx = atom.GetIdx()

    smiles = Chem.MolToSmiles(mol_pdb,True,False,idx)
    print("Ligand Smiles: ",smiles)
    map = mol_pdb.GetProp('_smilesAtomOutputOrder')
    atomorder = [int(i) for i in map.replace('[','').replace(']','').split(',') if i != '''''']
    orderarray = np.array(atomorder)
    indx_array = np.empty(num_atoms)  # index array
    for i in range(num_atoms):
        indx = int(orderarray[i])
#        print(i,indx,mol_pdb.GetAtomWithIdx(indx).GetPDBResidueInfo().GetName())
        indx_array[i] = indx

    mol = Chem.MolFromSmiles(smiles)
    rot_atom_pairs = mol.GetSubstructMatches(RotatableBondSmarts)
    print(rot_atom_pairs)
    ntors = len(rot_atom_pairs)
    k = 0
    for i in range(len(rot_atom_pairs)):
        indx1 = int(indx_array[ rot_atom_pairs[i][0] ])
        indx2 = int(indx_array[ rot_atom_pairs[i][1] ])
        atom1 = mol_pdb.GetAtomWithIdx(indx1).GetPDBResidueInfo().GetName()
        atom2 = mol_pdb.GetAtomWithIdx(indx2).GetPDBResidueInfo().GetName()
        iatom1 = mol_pdb.GetAtomWithIdx(indx1)
        iatom2 = mol_pdb.GetAtomWithIdx(indx2)
        #print(atom1, atom2, iatom1, iatom2)
        conf=mol_pdb.GetConformer(0)
        # Check if any angle is linear (torsion not defined). Using 10 degrees as cutoff.
        tors_flag = True
        for nn in iatom1.GetNeighbors():
           indx3 = nn.GetIdx()
           if (indx3 != indx2):
              angle = rdMolTransforms.GetAngleDeg(conf,indx3, indx1, indx2)
              #print(i, indx1, indx2, indx3, angle)
              if (abs(180.0 - abs(angle)) < 10.0):
                 tors_flag = False
        for nn in iatom2.GetNeighbors():
           indx3 = nn.GetIdx()
           if (indx3 != indx1):
              angle = rdMolTransforms.GetAngleDeg(conf,indx1, indx2, indx3)
              #print(i, indx1, indx2, indx3, angle)
              if (abs(180.0 - abs(angle)) < 10.0):
                 tors_flag = False
        if (tors_flag):
#           j = k + 3
           j = k + ncov + 1 #RS
           line = "move add mvtp tors weight 1.00 dmax 180.0 label tr" + str(j) + " fewer 0 sele atom ligand_@currligand 1 " \
                   + str(atom1) + " show end -\n"
           file.write(line)
           line = "                                                            sele atom ligand_@currligand 1 " \
                  + str(atom2) + " show end\n"
           file.write(line)
           k += 1
        else:
           ntors -= 1
    line = "\nincr ntors by " + str(ntors) + "\n\n"
    file.write(line)
#    print_part2(file)
    # Link all torsional moves (if smartmc is set by user)
# Remove for CHARMM-GUI until all MC bugs cleared
#    line = "if smartmc .eq. true then\n"
#    for i in range(ncov+1, ntors+1):
#       tr1 = "tr" + str(i)
#       tr2 = "tr" + str(i+1)
#       line += "   move link lab1 " + tr1 +  " lab2 " + tr2 + "\n"
#    line += "endif\n\n"
#    file.write(line)

    print_part2(file)
    file.close()

    # Write the mc_delete.str file for EnzyDock

    file = open(mcfile2,"w")
    print("Writing file: ", mcfile2)

    print_part3(file)
    # Unlink torsional moves (in reverse order)
# Remove for CHARMM-GUI until all MC bugs cleared
#    line = "if smartmc .eq. true then\n"
#    for i in range(ntors, ncov, -1):
#       tr1 = "tr" + str(i)
#       line += "   move link lab1 " + tr1 + "\n"
#    line += "endif\n\n"
#    file.write(line)

    print_part4(file)

    for i in range(ntors):
#        j = i + 3
        j = i + ncov + 1 #RS
        line = "move dele label tr" + str(j) + "\n"
        file.write(line)

    line = "\ndecr ntors by " + str(ntors) + "\n"
    file.write(line)
    print_part5(file)
    file.close()

#    print(Chem.MolToMolBlock(mol_pdb))
#    print(Chem.MolToMolBlock(mol))

    map = mol_pdb.GetProp('_smilesAtomOutputOrder')
    atomorder = [int(i) for i in map.replace('[','').replace(']','').split(',') if i != '''''']
    orderarray = np.array(atomorder)
    for i in range(num_atoms):
        indx = int(orderarray[i])
#        print(i,indx,mol_pdb.GetAtomWithIdx(indx).GetPDBResidueInfo().GetName())

def main():
    write_header()
    smiles_tors()
    write_footer()

if __name__ == "__main__":
    main()


