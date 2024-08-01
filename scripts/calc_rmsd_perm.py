#!/home/charmm-gui/local/miniconda3/bin/python
# Original version written by DTM 2020
# Permutations included by RS 03/2021
# Copyright © 2022 Dan T. Major

import copy
import math
import sys
import rdkit
from rdkit import Chem


'''
In this module there is a set of functions to evaluate RMSD between two
aligned structures of a molecule. The code can take into account
rotations around bonds (and maybe also overall rotation) that replaces 
chemically equivalent heavy atoms.
I.e. : if a t-butyl group was rotated by 60 degrees, despite the different indices
the code can consider each methyl group as replaced with another one.
Strategy: get lists of equivalent atoms using rdkit function "GetSubstructMatches".
Filter out the permutations that doesn't involve a mirror image (like two geminal positions,
or the two oxygen atoms in this molecule CH3-SO2-CH2CH3) by checking that the bond is rotatable
and the hybridization: 3 equivalent groups on SP3 atom, and 2 on SP2 atom.

Warning:
   The code excpects the first pdb file to have a CONECT section with multiple bonds
   written as duplicate CONECT records. No permutations can be defined with rdkit without
   the proper CONECT section.

Usage:
    rmsd_permutations_simplified.py pdb1 pdb2
For verbose mode (print used permutations):
    rmsd_permutations_simplified.py pdb1 pdb2 verbose
'''


### --- From here: from Dan's simplified rmsd script ---- ###


sys.path.append(".")

def write_header():
    print("\nRMSD program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

def write_footer():
    print("\nEnd RMSD program, returning to EnzyDock main\n")

def handle_input():
    # Handle input arguments and return values
    input_arg = sys.argv
    arg1 = str(input_arg[1])       # pdb file list
    arg2 = str(input_arg[2])       # pdb file list
    arg3 = False
    if len(sys.argv) > 3:
        arg = str(input_arg[3]).lower()
        if (arg == "verbose" or arg == "-verbose" or arg == "v" or arg == "-v"):
            arg3 = True
    return arg1, arg2, arg3


### --- From here: permutations auxiliary functions part ---- ###


def l_counter(i_list):
    # local counter: create a dictionary in which each key has a value 
    # according to the number of occurrences of the key in the list
    # written to avoid "from collections import Counter"
    d_counter = {}
    for item in i_list:
        if item in d_counter.keys():
            d_counter[item] += 1
        else:
            d_counter[item] = 1
    return d_counter

def find_exclude_bonded_inds(i_atom, molecule, exclude):
    # return a list of bounded atoms that do not appear in the exclude list
    bonded_inds = []
    bonded_atoms = molecule.GetAtoms()[i_atom].GetNeighbors()
    for i in range(len(bonded_atoms)):
        idx = bonded_atoms[i].GetIdx()
        if not (idx in exclude):
            bonded_inds.append(idx)
    return bonded_inds


### --- From here: permutations CHECK functions part  ---- ###


def is_bond_approved(atom_inx, permuted_atoms, molecule):
    # verify that bond order around permutable atoms is less/equal 1
    not_permutable_neighbors_inds = find_exclude_bonded_inds(atom_inx, molecule, permuted_atoms)
    if len(not_permutable_neighbors_inds) == 0:
        # meaning - a central atom (a terminal atom should be permutable)
        return True
    check = len(not_permutable_neighbors_inds)
    for neighbor in not_permutable_neighbors_inds:
        tested_bond = molecule.GetBondBetweenAtoms(atom_inx, neighbor)
        # bond order more than 1 - one cannot rotate around this
        ### print(tested_bond.GetBondTypeAsDouble())
        if tested_bond.GetBondTypeAsDouble() > 1.0:
            # what about the others (if there)
            check -= 1
    if check > 0:
        return True
    return False

def is_atom_approved(atom_inx, occurrences, molecule):
    # check for central atom hybridization and number of permuted atoms
    hybridization = str(molecule.GetAtoms()[atom_inx].GetHybridization())
    if hybridization == 'OTHER' or hybridization == 'UNSPECIFIED':
        #print("we have a problem: unrecognized hybridization")
        return False
    elif hybridization == 'S':
        #print("is it make sense that atom index " + str(atom_inx) + "has S hybridization?")
        return False
    elif hybridization == 'SP' and occurrences == 2:
        return True
    elif hybridization == 'SP2' and occurrences == 2:
        return True
    elif hybridization == 'SP3' and occurrences >= 3:
        return True
    # umbrella inversion around nitrogen 
    elif molecule.GetAtoms()[atom_inx].GetAtomicNum() == 7 and hybridization == 'SP3' and occurrences == 2:
        # avoid ammonium cases
        if  molecule.GetAtoms()[atom_inx].GetDegree() < 4:
            return True
    elif hybridization == 'SP3D' or hybridization == 'SP3D2':
        #print("too complicated case for atom " + str (atom_inx) +
        #     "with " + str(hybridization) + " hybridization. permutations around ignored")
        return False
    return False

def check_permutation_centeral_atom(neighbor, permuted_atoms, molecule):
    check_bond = is_bond_approved(neighbor[0], permuted_atoms, molecule)
    check_hybr = is_atom_approved(*neighbor, molecule)
    return check_bond and check_hybr


### --- From here: permutations main functions part ---- ###


def find_primary_permutations(mol):
    #identify groups of atoms that are replaced together
    
    permutations = mol.GetSubstructMatches(mol,uniquify=False)
    short_permutations_dicts = []
    # find out permutations of atoms (ignore atoms that are not permutated) - keep as dicts
    for permutation_i in permutations:
        short_d = {}
        for i, item in enumerate(permutation_i):
            if i != item:
                short_d[i] = item
        if short_d != {}:
            short_permutations_dicts.append(copy.deepcopy(short_d))
    # convert them into sets to easily filter out the primary replacements
    short_permutations_sets = []
    for short_d in short_permutations_dicts:
        short_permutations_sets.append(copy.deepcopy(set([(k, v) for k, v in short_d.items()])))
    # identify primary replacements
    primary = [True] * len(short_permutations_dicts)
    for i in range(len(short_permutations_sets)):
        for j in range(i + 1, len(short_permutations_sets)):
            if short_permutations_sets[i].issubset(short_permutations_sets[j]):
                primary[j] = False
            elif short_permutations_sets[j].issubset(short_permutations_sets[i]):
                primary[i] = False
    # create a list of the dictionaries describing the primary replacements
    primary_permutations_dict = []
    for i, short_d in enumerate(short_permutations_dicts):
        if primary[i]:
            primary_permutations_dict.append(copy.deepcopy(short_d))
    return primary_permutations_dict


def find_approved_permutations(molecule, primary_permutations_dict):
    '''
    find the atom that all 2 or 3 atoms are connected to and check its hybridization.
     accept the hybridization if it's sp2 that is connected to 2 permutable atoms,
     or if it's sp3 that is connected to 3 permutable atoms.
     REJECT if it's sp3 that is connected to 2 permutable atoms.
     also check for geminal position (that bond order of rotated bond is not greater than 1)
    '''

    # create list of the atoms bound to each permutable atom
    # later the list would be sorted to identify atoms that link between permutable atoms
    approved = [True] * len(primary_permutations_dict)
    neighbors = []
    for j, primary_permutation in enumerate(primary_permutations_dict):
        bonded_inds = []
        for i_atom in primary_permutation.keys():
            #do not include atoms that are permuted (in this permutation) to the list
            bonded_inds += find_exclude_bonded_inds(i_atom, molecule, primary_permutation.keys())
        neighbors.append(copy.deepcopy(bonded_inds))
    for i, primary_permutation in enumerate(neighbors):
        # counts occurrence of atoms in the neighbors
        permutation_neighbors_count = l_counter(primary_permutation)
        # check which neighbor occurs more than once and check for permutation around it
        for neighbor in permutation_neighbors_count.items():
            if neighbor[1] > 1: # central atom
                check = check_permutation_centeral_atom(neighbor, primary_permutations_dict[i].keys(), molecule)
                approved[i] = approved[i] and check
    #print(approved)
    approved_permutations_dict = []
    for i, primary_permutation in enumerate(primary_permutations_dict):
        if approved[i]:
            approved_permutations_dict.append(copy.deepcopy(primary_permutation))
    return approved_permutations_dict, approved


### --- From here: RMSD part ---- ###


def get_poses(pose1_pdb, pose2_pdb):
    # read pdb files and make sure they contain the same molecule
    poses = []
    mol_tmp1 = Chem.MolFromPDBFile(pose1_pdb)
    mol_tmp2 = Chem.MolFromPDBFile(pose2_pdb)
    if ((mol_tmp1 != None) and (mol_tmp2 != None)):
       poses += [ Chem.MolFromPDBFile(pose1_pdb) ]
       poses += [ Chem.MolFromPDBFile(pose2_pdb) ]
    else:
       print("Warning: Problems with molecules %s %s ... " % (pose1_pdb, pose2_pdb))
    #validation
    if not(Chem.MolToSmiles(poses[0]) == Chem.MolToSmiles(poses[1])):
        print("Warning: Poses of different molecules:\n", Chem.MolToSmiles(poses[0]), "\n",
              Chem.MolToSmiles(poses[1]))
        #      Chem.MolToSmiles(poses[1]), "\n. exit")
        #exit()
    return poses

def merge_permutations(old_d, new_d):
    # merge two dictionaries. When conflict, use the second one
    merged_d = copy.deepcopy(old_d)
    for key in new_d.keys():
        merged_d[key] = new_d[key]
    return merged_d

def opt_rmsd(mol1, mol2, approved_permutations_dict, verbose):
    # get best rmsd with optional permutations
    min_rmsd = get_rmsd_dict(mol1, mol2, {}) # initial values without any permutations
    all_permutations = {}
    for permutation in approved_permutations_dict:
        new_permutations = merge_permutations(all_permutations, permutation)
        new_rmsd = get_rmsd_dict(mol1, mol2, new_permutations)
        if  new_rmsd < min_rmsd:
            min_rmsd = new_rmsd
            all_permutations = copy.deepcopy(new_permutations)
    if verbose:
        if approved_permutations_dict == {}:
            print("No available permutations for this molecule")
        elif all_permutations == {}:
            print("Permutations were checked but weren't needed in this case")
        else:
            print("All permutations used:\n", all_permutations)
    else:
        if approved_permutations_dict == {}:
            print("No available permutations for this molecule")
        elif all_permutations == {}:
            print("Permutations were checked but weren't needed in this case")
        else:
            print("Some permutations were used in this case")
    return min_rmsd

def get_rmsd_dict(mol1, mol2, permutation):
    # Calculate RMSD between two poses using a single permutation
    # (send in an empty dictionary for no permutation)
    atom_count1 = mol1.GetNumAtoms()
    sum_2 = 0.0
    for k in range(atom_count1):
        if k in permutation.keys():
            j = permutation[k]
        else:
            j = k
        pos1 = mol1.GetConformer().GetAtomPosition(k)
        pos2 = mol2.GetConformer().GetAtomPosition(j)
        sum_2 += (pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2
    rmsd = math.sqrt(sum_2 / atom_count1)
    #print("RMSD: %.3f" % (rmsd))
    return rmsd


### --- From here: "main" part ---- ###


def main():
    write_header()
    pose1_pdb, pose2_pdb, verbose = handle_input()
    #read in structures
    poses = get_poses(pose1_pdb, pose2_pdb)
    molecule = poses[0]
    # filter out primary permutations, both mirror and axial
    primary_permutations_dict = find_primary_permutations(molecule)
    #filter out approved (probably axial) permutations
    approved_permutations_dict, approved = find_approved_permutations(molecule, primary_permutations_dict)   
    #print("all available primary permutations:\n", primary_permutations_dict)
    #print("\nThey were approved (True) or rejected (False) as follows:\n", approved)
    #print("\nOnly approved permutations are (see cartoon above for details):\n", approved_permutations_dict)
    rmsd_no_permutation = opt_rmsd(poses[0], poses[1], {}, verbose)
    print("rmsd without permutations: %.3f" % (rmsd_no_permutation))
    rmsd_with_permutation = opt_rmsd(poses[0], poses[1], approved_permutations_dict, verbose)
    print("rmsd with    permutations: %.3f" % (rmsd_with_permutation))
    write_footer()

if __name__ == "__main__":
    main()
    



