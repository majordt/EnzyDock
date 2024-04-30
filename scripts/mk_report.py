#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9

import sys
import os
import numpy as np

def handle_input():
    # Handle input arguments and return values
    # set 0 @numligands
    # set 1 @runid
    # set 2 @cluster
    input_arg = sys.argv
    resdir = "../results/"
    numligands = int(input_arg[1])
    runid = input_arg[2].lower() # CHARMM sends in upper case parameters, change to lower
    cluster = input_arg[3].lower()
    
    return resdir, numligands, runid, cluster

def parse_efile(fname, iskip):
    with open(fname, 'r') as f:
        lines = f.readlines()
        l1, l2, l3, l4 = [],[],[],[]
        if len(lines) <= iskip: # iskip is the number of lines in the header
            return l1, l2, l3, l4 # empty lists
        for line in lines[iskip:]:
            e1, e2, e3, cl = line.split()
            l1.append(float(e1))
            l2.append(float(e2))
            l3.append(float(e3))
            l4.append(cl)
    return l1, l2, l3, l4


def score(e_lists, factors):
    if (len(e_lists) == len(factors)):
        if (len (e_lists) > 1):
            e_lists_np = np.array(e_lists)
            factors_np = np.array(factors)
            iscore = np.matmul(e_lists_np.T, factors_np)
            return iscore.tolist()
        else:
            iscore = []
            for i in e_lists[0]:
                iscore.append(i*factors[0])
            return (iscore)
    else:
        print("number of factors specified doesn't match number of lists, assuming similar weight to all")
        factors = []
        for i in e_lists:
            factors.append(1)
        e_lists_np = np.array(e_lists)
        factors_np = np.array(factors)
        iscore = np.matmul(e_lists_np.T, factors_np)
        return iscore.tolist()

def write_report(iscore, clust_mm, resdir, currligand, runid):
    fname = resdir + str(currligand) + "_final_report_" + runid + ".dat"
    with open(fname, 'w') as f:
        f.write("*\n  SUMMARY FOR LIGAND NUMBER "+str(currligand)+"\n")
        f.write("ISCORE    CLUSTMITER\n")
        f.write("====================\n")
        for s, cl in zip(iscore, clust_mm):
            f.write(str(round(s,3)) + "    " + cl + "\n")

def write_header():
    print("\nCreating a report for ranking poses\n")

def write_footer():
    print("\nFininshed creating a report for ranking poses, returning to EnzyDock main\n")

def main():
    write_header()
    resdir, numligands, runid, cluster = handle_input()
    for currligand in range(1,numligands+1):
        
        # parse files
        fname = resdir + str(currligand) + "_summary_mm_" + runid + ".dat"
        mm = False
        if os.path.exists(fname):
            # etot: total ligand energy
            # intere: interaction energy + restraint energy
            # intere0: interaction energy
            etot_mm,intere_mm,intere0_mm,clust_mm = parse_efile(fname, 10)
            if len(etot_mm) > 0:
                mm = True
                #print("etot_mm,intere_mm,intere0_mm,clust_mm")
                #print(etot_mm,intere_mm,intere0_mm,clust_mm,sep="\n")
            
        
        # if gbsw is called with explicit water it is ignored
        fname = resdir + str(currligand) + "_c_summary_gbsw_" + runid + ".dat"
        gbsw = False
        if os.path.exists(fname):
            # etot (complex): Total GBSW energy (includes all terms)
            # elstat (complex): GBSW electrostatic energy
            # vdw (complex): GBSW vdW energy
            etot_gb,elstat_gb,vdw_gb,clust_gb = parse_efile(fname, 9)
            if len(etot_gb) > 0:
                gbsw = True
                #print("etot_gb,elstat_gb,vdw_gb,clust_gb")
                #print(etot_gb,elstat_gb,vdw_gb,clust_gb,sep="\n")
        
        # if pbeq is called with explicit water it is ignored
        fname = resdir + str(currligand) + "_c_summary_pbeq_" + runid + ".dat"
        pbeq = False
        if os.path.exists(fname):
            # dener: Electrostatic solvation energy difference (P+L-->PL)
            # dener80: PL interaction energy in water
            # dener1: PL interaction energy in vacuum
            dener_pb,dener80_pb,dener1_pb,clust_pb = parse_efile(fname, 8)
            if len(dener_pb) > 0:
                pbeq = True
                #print("dener_pb,dener80_pb,dener1_pb,clust_pb")
                #print(dener_pb,dener80_pb,dener1_pb,clust_pb,sep="\n")
       
        # if qm/mm is called with different number of atoms it is ignored
        fname = resdir + str(currligand) + "_summary_qmmm_" + runid + ".dat"
        qmmm = False
        if os.path.exists(fname):
            # qm/mm(tot): qm/mm potential energy
            # qm/mm(elstat): electrostatic part of qm/mm energy
            # qm/mm(vdw): vdw part of qm/mm energy
            etot_qm,elstat_qm,vdw_qm,clust_qm = parse_efile(fname, 10)
            if len(etot_qm) > 0:
                qmmm = True
                #print("etot_qm,elstat_qm,vdw_qm,clust_qm")
                #print(etot_qm,elstat_qm,vdw_qm,clust_qm,sep="\n")
        
        # integrate energies (needs refinement + testing)
        # gbsw pbeq and qmmm are not calculated for covalent ligands (@lastinter)
        # so in one system some ligands might have these energies and others not
        if gbsw == True:    
            e_lists = [etot_mm, elstat_gb, vdw_gb]
            factors = [1.0, 0.1, 1.0]

        else:
            e_lists = [etot_mm]
            factors = [1.0]
        iscore = score(e_lists, factors)
        write_report(iscore, clust_mm, resdir, currligand, runid)
    write_footer()
     
if __name__ == "__main__":
    main()
