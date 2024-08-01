#!/home/charmm-gui/local/miniconda3/bin/python
# RS: superpose_cluster_ligands as part of EnzyDock
# RS: compare to cluster_ligands.py, aviod duplicates of functions. Include symmetry.
# To test as standalone use:
# superpose_cluster_ligands.py 1 abyu 100 1.0
# Copyright © 2022 Dan T. Major

import math
import sys
import os
from read_write_pdb import Atom, PDB
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina
# RS
from rdkit.Chem import TorsionFingerprints
# avoid duplicates
from cluster_ligands import write_file, get_rmsd, cluster, write_footer, raise_nofile_err

def write_header():
    print("\nLigand clustering program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

def handle_input():
    # Handle input arguments and return values
    # set 0 @currligand
    # set 1 @runid
    # set 2 @maxconf
    # set 3 @clusterwidth

    input_arg = sys.argv
    currligstr = input_arg[1].lower() # CHARMM sends in upper case parameters, change to lower
    runid = input_arg[2].lower() # CHARMM sends in upper case parameters, change to lower
    bind = "mdmc_min"
    pre = "0"
    maxconf = int(input_arg[3])
    distTresh = float(input_arg[4])
    # return arg0, arg1, arg2, maxit, maxrot, distTresh, numseg
    return currligstr, runid, bind, pre, maxconf, distTresh

# def read_pdb_files(str0, str1, iligand, inum, n):
def read_pdb_files(currligstr, runid, bind, pre, maxconf):
    # Read pdb files and return list of molecular objects, directories and files
    print("Opening pdb files...\n")
    results_dir = "../pdb"
    # errfile = open('error.log', 'w')
    # try:
        # os.mkdir(results_dir)
    # except OSError as error:
        # print(error, file = errfile)
    mol = []
    file_list = []
    dirfile_list = []
    pdb_file0 = currligstr + "_" + runid + "_" + bind
    skip = 0
    for conf in range(maxconf):
         # @currligand_@runid_mdmc_min_@conf_0.crd
         pdb_file = pdb_file0 + "_" + str(conf+1) + "_" + pre + ".pdb"
         print(pdb_file)
         #file_list += [ pdb_file ]
         pdb_file = results_dir + "/" + pdb_file
         #dirfile_list += [ pdb_file ]
         subprocess.call(["/bin/rm","-f","temp"])
         tmpfilename = 'temp'
         pdb = PDB(pdb_file, tmpfilename)
         atoms = pdb.get_atoms(to_dict=False)
         pdb.add_element()
         pdb.write_pdb(tmpfilename)
         subprocess.call(["/bin/mv","temp",pdb_file])
         #mol += [ Chem.MolFromPDBFile(pdb_file) ]
         mol_tmp = Chem.MolFromPDBFile(pdb_file)
         if mol_tmp != None:
            mol += [ mol_tmp ]
            # @currligand_@runid_mdmc_min_@conf_0.crd
            pdb_file = pdb_file0 + "_" + str(conf+1) + "_" + pre + ".pdb"
            file_list += [ pdb_file ]
            pdb_file = results_dir + "/" + pdb_file
            dirfile_list += [ pdb_file ]
         else:
            print("Problems with molecule %s ... skipping to next" % pdb_file)
            skip += 1
    if skip > 0:
       print("Skipped %d out of %d structures" % (skip, maxconf))
       if skip == maxconf:
          raise_nofile_err(currligstr)
    print("\n")
    return mol, results_dir, file_list, dirfile_list

# def get_rmsd(mol, align=False):
def get_tfd(mol):
    # Calculate TFD between all simulated conformers and return triangular matrix of values
    # in format required by clustering function
    mol_count = len(mol)
    print("\nNumber of molecules: %d\n" % mol_count)

    dm = []
    # Calculate TFD for all ligand states
    # For loop details, see: https://www.rdkit.org/docs/source/rdkit.ML.Cluster.Butina.html
    for i in range(mol_count):
        for j in range(i):
            # calculate TFD and store in an array
            tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol[i], mol[j])
            # tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol[i].GetConformer(), mol[j].GetConformer())
            # atom_count = mol[i].GetNumAtoms()
            # sum_2 = 0.0
            # for k in range(atom_count):
                # pos1 = mol[i].GetConformer().GetAtomPosition(k)
                # pos2 = mol[j].GetConformer().GetAtomPosition(k)
                # sum_2 += (pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2
            # rmsd = math.sqrt(sum_2 / atom_count)
            print("Comparing molecule #%4d with molecule #%4d (%4d molecules in all). TFD: %.3f" % (i, j, mol_count,tfd))
            #print("RMSD: %.3f" % rmsd)
            dm.append(tfd)
    return dm, mol_count

def file_handle(cs, results_dir, file_list, dirfile_list, currligstr, runid, bind, pre):
    # Copy files into cluster directories and identify lowest energy pose per cluster
    # This function makes sure that all files are in place for QM calculations
    # on the lowest energy representative of each cluster

    debug = False
    errfile = open('error.log', 'w') 
    lig_dir = "ligand_" + currligstr
    # Clean up old cluster directories first
    lig_results_dir = results_dir + "/" + lig_dir
    cluster_dir = lig_results_dir + "/" + runid + "_" + bind + "_cluster" + "*" + "/"
    command = "rm -rf " + cluster_dir
    print("\nCleaning up old cluster directories: ", command, "\n")
    os.system(command)
    print("PDB files clustered:\n")
    try:
       os.mkdir(lig_results_dir)
    except OSError as error:
       print(error, file = errfile) 
    for i in range(len(cs)):
        cluster_dir = results_dir + "/" + lig_dir + "/" + runid + "_" + bind + "_cluster" + str(i+1) + "/"
        try:
           os.mkdir(cluster_dir)
        except OSError as error:
           print(error, file = errfile)
        min_energy = 1.0e20   # Large number
        for j in range(len(cs[i])):
            k = cs[i][j]
            print("Cluster: %4d, %4d" % (i+1, k), " file name: %s" % dirfile_list[k])
            # Copy ligand poses into correct cluster directory
            command = "cp -f " + dirfile_list[k] + " " + cluster_dir
            os.system(command)
            file = open(dirfile_list[k], "r")
            line = file.readline()
            energy = float(line.split(' ')[5])   # This is the energy term
            if (energy < min_energy):
                min_energy = energy
                min_file = file_list[k]
        # First cp copies minimum ligand pose of cluster excluding flexible residues
        #command = "cp -f " + results_dir + "/" + min_file + " " + cluster_dir + "min_clust" + str(i+1) + "_" + min_file
        command = "cp -f " + results_dir + "/" + min_file + " " + cluster_dir + "ligand_" + currligstr + "_min_clust" + str(i+1) + ".pdb"
        try:
            os.system(command)
        except OSError as error:
            print(error, file = errfile)
        # copy + rename minimum energy structure to match dock_ligand.inp expectations
        # pdb file:
        # @resDIR/@currligand_@runid_@iter_0.pdb
        targetn = currligstr + "_" + runid + "_" + str(i+1) + "_" + pre + ".crd"
        targetdir = "../results"
        command = "cp -f " + results_dir + "/" + min_file + " " + targetdir + "/" + targetn
        try:
            os.system(command)
        except OSError as error:
            print(error, file = errfile)
        # crd file:
        # @resDIR/@currligand_@runid_@iter_0.crd
        min_file = min_file.rstrip("pdb") + "crd"
        targetn = currligstr + "_" + runid + "_" + str(i+1) + "_" + pre + ".crd"
        targetdir = "../results"
        command = "cp -f " + results_dir + "/" + min_file + " " + targetdir + "/" + targetn
        try:
            os.system(command)
        except OSError as error:
            print(error, file = errfile)
        for j in range(len(cs[i])):
            k = cs[i][j]
            # remove ligand poses into correct cluster directory
            command = "rm -f " + dirfile_list[k]
            os.system(command)
            # remove ligand poses into correct cluster directory
            crdname=dirfile_list[k][:-3]+'crd'
            command = "mv " + crdname + " " + cluster_dir
            os.system(command)
    
def tar_all(currligstr, bind):
    results_dir = "../pdb"
    # tar -zcf --remove-files file.tar.gz /path/to/dir/
    #command = "tar -zcf --remove-files "+results_dir+"/lig"+currligstr+".tar.gz "+results_dir+"/*"+bind+"* "+results_dir+"/ligand_"+currligstr
    os.system("pwd")
    idir=os.getcwd()
    idir=idir+"/"+results_dir
    command = "tar -zc --remove-files -C "+idir+" --file="+idir+"/lig"+currligstr+".tar.gz ligand_"+currligstr
    print(command,end="\n\n")
    try:
        os.system(command)
    except OSError as error:
        print(error, file = errfile)
    os.system("pwd")

def main():
    write_header()
    # str0, str1, iligand, inum, distTresh = handle_input()
    currligstr, runid, bind, pre, maxconf, distTresh = handle_input()
    mol, results_dir, file_list, dirfile_list = read_pdb_files(currligstr, runid, bind, pre, maxconf)
    # dm, mol_count = get_rmsd(mol, align=True)
    dm, mol_count = get_tfd(mol)
    cs = cluster(dm, mol_count, distTresh, currligstr, fname='niters', charmmvar='maxit')
    # file_handle(cs, results_dir, file_list, dirfile_list, str0, str1, n)
    file_handle(cs, results_dir, file_list, dirfile_list, currligstr, runid, bind, pre)
    tar_all(currligstr, bind)
    write_footer()

if __name__ == "__main__":
    main()
    
