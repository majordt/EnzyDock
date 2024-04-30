#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
# To test as standalone use:
# cluster_ligands_perm.py 1 cotb2_wt 5 5 1.0 1
# 
# Warnings:
#   - For permutations to be defined in rdkit, input pdb file (@pdbDIR/ligand_@currligand.pdb)
#     must contain CONECT section with mulitple bonds written as duplicate CONECT records.
#   - Permutations are blind to patches done after reading the input pdb file mentioned above.
#     I.e.: if a cofactor or the protein is connected to a symmetric position, this will (probably)
#     break the symmetry BUT the code cannot detect this.
# Original version written by DTM 2020
# Permutations included by RS 03/2021
# Copyright © 2022 Dan T. Major

import math
import sys
import os
sys.path.append(".")
from read_write_pdb import Atom, PDB
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina

#RS:
import calc_rmsd_perm as per_rmsd
# avoid duplicates
from cluster_ligands import write_header, write_file, handle_input, raise_nofile_err
from cluster_ligands import read_pdb_files, cluster, file_handle, write_footer
#RS


#def write_header():
#    print("\nLigand clustering program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

# def write_file(file):
    # file.write(
# '''* Number of clusters
# *
# set ncluster = ''')

# def handle_input():
    # # Handle input arguments and return values
    # input_arg = sys.argv
    # arg0 = "lig_"
    # arg1 = input_arg[1]
    # arg2 = input_arg[2].lower()   # CHARMM sends in upper case parameters, change to lower
    # maxit = int(input_arg[3])
    # maxrot = int(input_arg[4])
    # distTresh = float(input_arg[5])
    # numseg = int(input_arg[6])
    # return arg0, arg1, arg2, maxit, maxrot, distTresh, numseg

# def read_pdb_files(arg0, arg1, arg2, maxit, maxrot):
    # # Read pdb files and return list of molecular objects, directories and files
    # print("Opening pdb files...\n")
    # results_dir = "../results/"
    # pdb_file0 = arg0 + arg1 + "_" + arg2 + "_"
    # mol = []
    # file_list = []
    # dirfile_list = []
    # for i in range(maxit):
        # for j in range (maxrot):
            # pdb_file = pdb_file0 + str(i+1) + "_" + str(j+1) + ".pdb"
            # print(pdb_file)
            # file_list += [ pdb_file ]
            # pdb_file = results_dir + pdb_file
            # dirfile_list += [ pdb_file ]
            # subprocess.call(["/bin/rm","-f","tmp"])
            # tmpfilename = 'tmp'
            # pdb = PDB(pdb_file, tmpfilename)
            # atoms = pdb.get_atoms(to_dict=False)
            # pdb.add_element()
            # pdb.write_pdb(tmpfilename)
            # subprocess.call(["/bin/mv","tmp",pdb_file])
            # # Hardwire with use of obabel to overcome problems with rdkit not recognizing
            # # element in 'atom name' field and CHARMM not printing element in element symbol field.
            # # DTM 25/01/2021 The need for obabel has been removed by independent code fixing PDB format for RDKit
            # #subprocess.call(["/bin/rm","-f","tmp"])
            # #subprocess.call(["obabel","-ipdb",pdb_file,"-opdb","tmp"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # #subprocess.call(["/bin/mv","tmp",pdb_file])
            # mol += [ Chem.MolFromPDBFile(pdb_file) ]
    # print("\n")
    # return mol, results_dir, file_list, dirfile_list

#RS:
def permutations(currligand):
    # read the molecule from a pdb file containing connect section to define permutations 
    connect_pdb_file = "../pdb/ligand_" + currligand + ".pdb"
    molecule = Chem.MolFromPDBFile(connect_pdb_file)
    # filter out primary permutations, both mirror and axial
    primary_permutations_dict = per_rmsd.find_primary_permutations(molecule)
    #filter out approved (probably axial) permutations
    approved_permutations_dict, approved = per_rmsd.find_approved_permutations(molecule, primary_permutations_dict)   
    return approved_permutations_dict

def get_rmsd(mol, approved_permutations_dict, use_permutation):
    # Calculate RMSD between all docked poses and return triangular matrix of values
    # in format required by clustering function
    mol_count = len(mol)
    print("\nNumber of molecules: %d\n" % mol_count)
    dm = []
    # Calculate RMSD for all ligand states (verified against rms_cur in Pymol)
    # For loop details, see: https://www.rdkit.org/docs/source/rdkit.ML.Cluster.Butina.html
    for i in range(mol_count):
        for j in range(i):
            #print("Aligning molecule #%d with molecule #%d (%d molecules in all)" % (i, j, mol_count))
            # calculate RMSD and store in an array
            if (use_permutation == False) or (approved_permutations_dict == {}):
                atom_count = mol[i].GetNumAtoms()
                sum_2 = 0.0
                for k in range(atom_count):
                    pos1 = mol[i].GetConformer().GetAtomPosition(k)
                    pos2 = mol[j].GetConformer().GetAtomPosition(k)
                    sum_2 += (pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2
                rmsd = math.sqrt(sum_2 / atom_count)
            else:
                rmsd = per_rmsd.opt_rmsd(mol[i], mol[j], approved_permutations_dict, False)
            print("Aligning molecule #%4d with molecule #%4d (%4d molecules in all). RMSD: %.3f" % (i, j, mol_count,rmsd))
            #print("RMSD: %.3f" % rmsd)
            dm.append(rmsd)
    return dm, mol_count
#RS

# def cluster(dm, mol_count, distTresh, arg1):
    # # Perform clustering of all ligands

    # #print(len(dm))
    # #print(dm)
    # cs = Butina.ClusterData(dm,mol_count,distTresh,isDistData=True,reordering=True)
    # print("\nNumber of clusters: %d\n" % len(cs))
    # # Write number of clusters to file for CHARMM to pick up
    # cfilename = 'nclusters' + arg1 + '.str'
    # nclusterfile = open(cfilename, 'w')
    # write_file(nclusterfile)
    # print(len(cs), file = nclusterfile)
    # k = 1
    # for c in cs:
        # print("Cluster %4d: " % k, '{}'.format(c))
        # k += 1
    # return cs

# def file_handle(cs, results_dir, file_list, dirfile_list, arg0, arg1, arg2, numseg):
    # # Copy files into cluster directories and identify lowest energy pose per cluster
    # # This function makes sure that all files are in place for MM and QM/MM calculations
    # # on the lowest energy representative of each cluster

    # debug = False
    # errfile = open('error.log', 'w') 
    # # Clean up old cluster directories first
    # cluster_dir = results_dir + arg1 + "_" + arg2 + "_" + "cluster"
    # command = "rm -rf " + cluster_dir + "*"
    # print("\nCleaning up old cluster directories: ", command, "\n")
    # os.system(command)
    # print("PDB files clustered:\n")
    # for i in range(len(cs)):
        # cluster_dir = results_dir + arg1 + "_" + arg2 + "_" + "cluster" + str(i+1) + "/"
        # try:
            # os.mkdir(cluster_dir)
        # except OSError as error:
            # print(error, file = errfile)
        # min_energy = 1.0e20   # Large number
        # for j in range(len(cs[i])):
            # k = cs[i][j]
            # print("Cluster: %4d, %4d" % (i+1, k), " file name: %s" % dirfile_list[k])
            # # Copy ligand poses into correct cluster directory
            # command = "cp -f " + dirfile_list[k] + " " + cluster_dir
            # os.system(command)
            # file = open(dirfile_list[k], "r")
            # line = file.readline()
            # energy = float(line.split(' ')[2])   # This is the 'inter' term, which includes restraints if any
            # if (energy < min_energy):
                # min_energy = energy
                # min_file = file_list[k]
        # # First cp copies minimum ligand pose of cluster excluding flexible residues
        # command = "cp -f " + results_dir + min_file + " " + cluster_dir + "min_clust" + str(i+1) + "_" + min_file
        # try:
            # os.system(command)
        # except OSError as error:
            # print(error, file = errfile)
        # # Second cp copies minimum ligand pose of cluster including flexible residues
        # if arg0 in command:
            # command = command.replace(arg0,"")
        # try:
            # os.system(command)
        # except OSError as error:
            # print(error, file = errfile)
        # # Copy crd file which will be used by MM and QM/MM rescoring calculations
# #    if arg1 in command:
        # command = "cp -f " + results_dir + min_file  + " " + cluster_dir + "min_clust" + str(i+1) + "_" + arg1 + "_" + arg2 + ".crd"
        # command = command.replace("pdb","crd")
# #    if arg0 in command:
# #        command = command.replace(arg0,"")
        # try:
            # os.system(command)
        # except OSError as error:
            # print(error, file = errfile)
        # # Copy flexible residues files, must be done per segment
        # for k in range(numseg):
            # in_flex_file = "flex" + str(k+1) + "_" + min_file
            # in_flex_file = in_flex_file.replace(arg0,"")
            # out_flex_file = "flex" + str(k+1) + "_min_clust" + str(i+1) + "_" + arg1 + "_" + arg2 + ".pdb"
            # in_flex_file = results_dir + in_flex_file
            # out_flex_file = cluster_dir + out_flex_file
            # command = "cp -f " + in_flex_file  + " " + out_flex_file
            # if os.path.isfile(in_flex_file):
                # os.system(command)
            # else:
                # if (debug): print("File %s doesn't exist" % in_flex_file)
        # # Copy waters
        # in_flex_file = "waterin_" + min_file
        # in_flex_file = in_flex_file.replace(arg0,"")
        # out_flex_file = "waterin_min_clust" + str(i+1) + "_" + arg1 + "_" + arg2 + ".pdb"
        # in_flex_file = results_dir + in_flex_file
        # out_flex_file = cluster_dir + out_flex_file
        # command = "cp -f " + in_flex_file  + " " + out_flex_file
        # if os.path.isfile(in_flex_file):
            # os.system(command)
        # else:
            # if (debug): print("File %s doesn't exist" % in_flex_file)

# def write_footer():
    # print("\nEnd clustering program, returning to EnzyDock\n")

def main():
    write_header()
    arg0, arg1, arg2, maxit, maxrot, distTresh, numseg = handle_input()
    mol, results_dir, file_list, dirfile_list = read_pdb_files(arg0, arg1, arg2, maxit, maxrot)
    #RS:
    approved_permutations_dict = permutations(arg1)
    dm, mol_count = get_rmsd(mol, approved_permutations_dict, use_permutation = True)
    #dm, mol_count = get_rmsd(mol)
    #RS
    cs = cluster(dm, mol_count, distTresh, arg1)
    file_handle(cs, results_dir, file_list, dirfile_list, arg0, arg1, arg2, numseg)
    write_footer()

if __name__ == "__main__":
    main()



