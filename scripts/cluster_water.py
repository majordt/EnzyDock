#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
# To test as standalone use:
# cluster_water.py <int> outfile.pdb
# 
# Clustering of water molecules using DBSCAN. Used by EnzyWater
# Original version written by DTM 08/2022
# Copyright © 2022 Dan T. Major

import math
import sys
import os
import random
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.neighbors import NearestNeighbors
from kneed import KneeLocator
sys.path.append(".")
from read_write_pdb import Atom, PDB
import subprocess
from rdkit import Chem

def write_header():
    print("\nWater clustering program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

def handle_input():
    ''' Handle input arguments and return values '''

    input_arg = sys.argv
    arg0 = "waterexpl_"
    arg1 = int(input_arg[1])
    arg2 = input_arg[2].lower()
    results_dir = "../pdb/"
    return arg0, arg1, arg2, results_dir

def read_pdb_files(arg0, results_dir, mcwat_ntot):
    ''' Read pdb files and return list of molecular objects, directories and files '''
    print("Opening pdb files...\n")
    pdb_file0 = arg0
    mol = []
    file_list = []
    dirfile_list = []
    rand_num = random.randint(0, sys.maxsize)
    tmpfilename = "tmp" + str(rand_num)
    skip = 0
    for i in range(mcwat_ntot):
         pdb_file = pdb_file0 + str(i+1) + ".pdb"
         print(pdb_file)
         #file_list += [ pdb_file ]
         pdb_file = results_dir + pdb_file
         #dirfile_list += [ pdb_file ]
         #subprocess.call(["/bin/rm","-f",tmpfilename])
         pdb = PDB(pdb_file, tmpfilename)
         atoms = pdb.get_atoms(to_dict=False)
         pdb.add_element()
         #pdb.write_pdb(tmpfilename)
         #subprocess.call(["/bin/mv",tmpfilename,pdb_file])
         # Catch problems with molecule (RDKit returns None if there's a problem)
         mol_tmp = Chem.MolFromPDBFile(pdb_file)
         if mol_tmp != None:
            mol += [ mol_tmp ]
            pdb_file = pdb_file0 + str(i+1) + ".pdb"
            file_list += [ pdb_file ]
            pdb_file = results_dir + pdb_file
            dirfile_list += [ pdb_file ]
         else:
            print("Problems with molecule %s ... skipping to next" % pdb_file)
            skip += 1
#         mol += [ Chem.MolFromPDBFile(pdb_file) ]
         if skip > 0:
            print("Skipped %d out of %d structures" % (skip, mcwat_ntot))
    return mol, results_dir, file_list, dirfile_list

def dbscan_cluster(mol, min_samples, eps):
    ''' Cluster water molecules using DBSCAN '''
#    print(len(mol))
    max_min_samples = len(mol)
    data = []
    num_waters = len(mol)
    for i in range(num_waters):
#       print(mol[i].GetNumAtoms())
       for j, atom in enumerate(mol[i].GetAtoms()):
          positions = mol[i].GetConformer().GetAtomPosition(j)
          #print(atom.GetSymbol(), positions.x, positions.y, positions.z)
          data += [[positions.x, positions.y, positions.z]]
    #print(data)
    # Try estimate epsilon using knee method
    nearest_neighbors = NearestNeighbors(n_neighbors=min_samples)
    neighbors = nearest_neighbors.fit(data)
    distances, indices = neighbors.kneighbors(data)
    distances = np.unique(distances[:,min_samples-1], axis=0) # Returns unique sorted list
    #for i in range(len(distances)):
    #   print(i,distances[i])
    i = np.arange(len(distances))
    knee = KneeLocator(i, distances, S=1.0, curve='convex', \
             direction='increasing', interp_method='polynomial')
    print("\nKnee found %d %f" % (knee.knee, distances[knee.knee]))
    # Hardwire sanity check on eps (shouldn't be too low or higher than normal H-bond distances)
    eps = distances[knee.knee]
    if eps > 2.8:
       eps = 2.8
       print("Knee too high, using upper eps value of %f" % eps)
    elif eps < 0.9:
       eps = 0.9
       print("Knee too low, using lower eps value of %f" % eps)
    eps = 0.9   # Optimized parameters
    min_samples = int(max_min_samples/2)   # Optimized parameters
    data = np.array(data)
    print("Clustering waters using DBSCAN with eps %f and min_samples %d" % (eps, min_samples))
    #print("GetAtoms", len(mol[0].GetAtoms()))
    model = DBSCAN(eps=eps, min_samples=min_samples)
    model.fit_predict(data)
    pred = model.fit_predict(data)
    labels = model.labels_

    #print("Number of clusters found: {}".format(len(set(labels))))
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    if (n_clusters_ > 1 and min_samples <= max_min_samples):
       print("Silhouette Coefficient: %0.3f"
             % metrics.silhouette_score(data, labels))
#       for i in range(len(model.labels_)):
#          print('Cluster for point: ', i, model.labels_[i])
    return(data, model, n_clusters_)

def file_handle(data, model, mol, results_dir, file_list, dirfile_list, arg0, outfile): 
    ''' Determine centroid of each cluster and write to file '''

    #debug = False
    #errfile = open('error.log', 'w') 
    # Clean up old cluster directories first
#    cluster_dir = results_dir + arg0 + "cluster"
#    command = "rm -rf " + cluster_dir + "*"
#    print("\nCleaning up old cluster directories: ", command, "\n")
#    os.system(command)
#    print(model)
    labels = model.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    # Get the frequency count of each cluster
    cluster_counts = np.bincount(labels[labels>=0])
#    print(cluster_counts)
#    print(np.shape(cluster_counts))
#    print(len(cluster_counts))
#    print(len(labels))
#    print(labels)
    # Compute clusters centroids
    sum_x = np.zeros(len(cluster_counts))
    sum_y = np.zeros(len(cluster_counts))
    sum_z = np.zeros(len(cluster_counts))
    for i in range(len(labels)):
        j = labels[i]
        if j >= 0:
           sum_x[j] += data[i][0]
           sum_y[j] += data[i][1]
           sum_z[j] += data[i][2]
#           print(data[i][0],data[i][1],data[i][2],sum_x[j],sum_y[j],sum_z[j])
#    print(np.shape(sum_x))
#    print(np.shape(cluster_counts))
    x_centroid = sum_x / cluster_counts
    y_centroid = sum_y / cluster_counts
    z_centroid = sum_z / cluster_counts
    # Compute occupation frequency for each water cluster
    sum_cluster_counts = sum(cluster_counts)
    cluster_counts = 100.0*cluster_counts/sum_cluster_counts
    # Print centroids and occupancies to PDB file
    tmpfilename = results_dir + outfile
    print("\nPDB file of clustered waters: %s" % tmpfilename)
    pdb = PDB(file2=tmpfilename)
    for i in range(len(cluster_counts)):
#        print(i, x_centroid[i], y_centroid[i], z_centroid[i], cluster_counts[i])
        pdb.add_atoms(serial=i+1, name='OH2', resName='TIP3', resSeq=i+1, x=x_centroid[i], \
                      y=y_centroid[i], z=z_centroid[i], occupancy=cluster_counts[i], element='O')
    pdb.write_pdb(tmpfilename)

def write_footer():
    print("\nEnd clustering program, returning to EnzyDock main\n")

def main():
    write_header()
    arg0, mcwat_ntot, outfile, results_dir = handle_input()
    mol, results_dir, file_list, dirfile_list = read_pdb_files(arg0, results_dir, mcwat_ntot)
    ## Iterate over DBSCAN until a reasonable cluster is found
    #n_clusters = 0
    #n_attempts = 1
    #max_attempts = 10
    #eps = 0.5
    #deps = 0.5
    min_samples = len(mol)   # Equal to number of water files
    #while (n_clusters <= 0 and n_attempts <= max_attempts):    
    data, model, n_clusters = dbscan_cluster(mol, min_samples, eps=None)
    #   eps += deps
    #   n_attempts += 1
    if (n_clusters > 0):
       file_handle(data, model, mol, results_dir, file_list, dirfile_list, arg0, outfile)
    else:
       print("\nEnzyWater WARNING: No clusters found...")
    write_footer()

if __name__ == "__main__":
    main()


