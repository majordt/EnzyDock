#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9

# This scripts enables the user to add noe restraints (not consensus) through a csv file.
# This script translates the csv to a CHARMM str file.
#
#Details:
#   * Atoms selected:
#        The user may choose any atom from the protein / cofactors / ligands.
#        Consensus selection (i.e. between two atoms from different ligands) CANNOT be done here.
#   * Implicit residues:
#        If the user chose implicit residues this script translates them to PNOE when the grid is on. 
#   * CSV File format:
#        The csv file needs to include the following columns:
#           segid_a,resid_a,type_a,segid_b,resid_b,type_b,rmin,rmax,k,ligands.
#        The "k" column is kmax. Fmax is set to twice its value. We tentatively suggest 25 or 50 as a default value for k.
#        The "ligands" column may contain a list (in python format, e.g: [1,3,5,4]) with the serial number of the ligands
#        for which the restrain should apply. If empty or "all" it is taken as all ligands.
#   * File location+name:
#         stream/consdef/userrestraints.csv. This script writes the stream/consdef/userrestraints.str file.
#   * Overwriting:
#        If a csv file in this name is found this script will OVERWRITE an existing userrestraints.str file.
#Usage:
#   proc_noes.py <numligands> <flexstr> <proteinname>
#   "flexstr" is a string that contain information about the flexible residues and the number of protein units.
#   It is created inside EnzyDock: @{proteinunit}_pep1_@{flex11}_@{flex12}_{...}_pep2_@{flex21}_@{flex22}_{...}

import subprocess
from os import path
import sys
import pandas as pd
from read_write_pdb import PDB #, Atom

def write_header():
    print("\nStrat: Create NOE/PNOE restraints from user's csv file\n")

def write_footer():
    print("\nEnd: Create NOE/PNOE restraints from user's csv file, returning to EnzyDock main\n")

def handle_input():
    # Handle input arguments and return values   
    # set 0 @numligands
    # set 1 @flexstr
    # set 2 @proteinname
    input_arg = sys.argv
    numligands = int(input_arg[1])
    flexstr = input_arg[2].lower()
    proteinname = input_arg[3].lower() # CHARMM sends uppercase
    return numligands, flexstr, proteinname

def process_flexstr(flexstr):
    flex_raw = flexstr.split("_")
    proteinunit = int(flex_raw[0])
    flex = []
    for i in flex_raw[1:]:
        if "pep" in i:
            flex.append([])
        else:
            flex[-1].append(int(i))
    # if len(flex) == proteinunit:
        # print("succeed")
    return flex, proteinunit

# process ligand condition
def process_lig_n(ligands, numligands):
    if pd.isna(ligands):
        ligands = ""
    if type(ligands) == str:
        if ligands.lower() in ["all",""]:
            ligands1 = []
        elif ligands.isnumeric():
            ligands1 = [int(ligands)]
        else:
            try:
                ligands1 = eval(ligands)
            except:
                line = "when processing stream/consdef/userrestraints.csv,\n"
                line += "   could not parse ligands definition, of "+str(ligands)+" takes it as \"all\""
                print(line)
                # write to enzydock.log file
                with open("enzydock.log",'a') as f:
                    f.write(line+"\n")
                ligands1=[]
            if set(ligands1) == set([i for i in range(1,numligands+1)]):
                ligands1=[]
            else:
                ligands1 = sorted(ligands1)
    elif type(ligands) == int:
        ligands1 = [ligands]
    else:
        line = "when processing stream/consdef/userrestraints.csv,\n"
        line += "   could not parse ligands definition, of "+str(ligands)+" takes it as \"all\""
        print(line)
        # write to enzydock.log file
        with open("enzydock.log",'a') as f:
            f.write(line+"\n")
        ligands1 = []
    return ligands1

def read_pdbf(proteinname,pep_unit):
    # better to read file once and pass atom list
    pdbdir = "../pdb/"
    pdb_file = pdbdir + proteinname + "_" + str(pep_unit) + ".pdb"
    print("reads",pdb_file, "...")
    pdb = PDB(pdb_file)
    atoms = pdb.get_atoms(to_dict=False)   
    return atoms

def read_csv(csvfname,numligands,proteinname,proteinunit):
    # idir="../stream/consdef/"
    # fname="userrestraints.csv"
    # df = pd.read_csv(idir + fname)
    df = pd.read_csv(csvfname)
    # lower case of segid and match ligand naming
    for seg in ["segid_a", "segid_b"]:
        df[seg] = df[seg].apply(lambda x: x.lower())
        df[seg] = df[seg].apply(lambda x: "ligand_@currligand" if "lig" in x.lower() else x)
    # take care of list of ligands for which noe should apply
    if ("ligands" not in df.columns) or (numligands == 1):
        df["ligands1"] = df.apply(lambda x: [], axis=1)
    else:
        df["ligands1"] = df.apply(lambda x: process_lig_n(x.ligands, numligands), axis=1)
    # read pdb files of protein if needed
    atoms = []
    for pep_unit in range(1,proteinunit+1):
        if pd.concat([df['segid_a'], df['segid_b']]).str.contains("pep"+str(pep_unit)).any():
            atoms.append(read_pdbf(proteinname,pep_unit))
        else:
            atoms.append("")
    return df, atoms

def stop_if_no_coor():
    with open("../stream/error_restraints.str", 'w') as f:
        lines = """**
*
stop
"""
        f.write(lines)

def get_coor(atoms,pep_unit,resid,itype):
    # better to read file once and pass atom list
    # pdbdir = "../pdb/"
    # pdb_file = pdbdir + proteinname + "_" + str(pep_unit) + ".pdb"
    # pdb = PDB(pdb_file)
    # atoms = pdb.get_atoms(to_dict=False)
    for atom in atoms[pep_unit-1]:
        if atom.resSeq == resid and atom.name.lower() == itype.lower():
            return atom.x, atom.y, atom.z
    line = "when processing stream/consdef/userrestraints.csv,\n"
    line += "   could not find atom pep"+str(pep_unit)+" "+str(resid)+" "+str(itype)+"\nTerminating EnzyDock run..."
    print(line)
    # write to enzydock.log file
    with open("enzydock.log",'a') as f:
        f.write(line+"\n")
    stop_if_no_coor()
    return None, None, None

def mk_condition(liglist):
    if len(liglist) == 0:
        line = ""
    else:
        # line = "if ("
        # for i in range(len(liglist)):
            # line += ("@currligand .eq. "+str(liglist[i]))
            # if i != len(liglist) - 1: # not the last ligand
                # line += " .or. "
                # if i%2 == 1: #add a new line after 2 conditions to avoid charmm long lines)
                    # line += "-\n"
        # line += ") then\n"
        line = "set donoe false\n"
        for i in range(len(liglist)):
            line += ("if @currligand .eq. "+str(liglist[i])+" set donoe true\n")
        line += "if @donoe .eq. true then\n"
    return line

def raise_sele_err(noe_idx):
    line = """
if ?nsel .eq. 0 then
   echo EnzyDock WARNING: bad atom definition at noe """ + str(noe_idx) + " for ligand_@currligand\n"
    line += """   ! Update readable log file
!   open append unit 33 form name enzydock.log
   write title unit 33
* EnzyDock WARNING: bad atom definition at noe """ + str(noe_idx) + " for ligand_@currligand\n"
    line += """* Terminating EnzyDock run...
*
   close unit 33
   stop
endif\n"""
    return line

def simple_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k):
    line = "noe\n   assign select atom "+str(segid_a)+" "+str(resid_a)+" "+str(type_a)+" show end -\n"
    line += "   select atom "+str(segid_b)+" "+str(resid_b)+" "+str(type_b)+" show end -\n"
    line += "   kmin 0.0 rmin "+str(rmin)+" kmax "+str(k)+" rmax "+str(rmax)+" fmax "+str(float(k)*2)+"\nend\n"
    return line

def simple_mk_pnoe(a_x, a_y, a_z, segid_b, resid_b, type_b, rmin, rmax, k):
    line = "noe\n   assign cnox "+str(a_x)+" cnoy "+str(a_y)+" cnoz "+str(a_z)+" -\n"
    line += "   select atom "+str(segid_b)+" "+str(resid_b)+" "+str(type_b)+" show end -\n"
    line += "   kmin 0.0 rmin "+str(rmin)+" kmax "+str(k)+" rmax "+str(rmax)+" fmax "+str(float(k)*2)+"\nend\n"
    return line

def pep_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k, flex, atoms):
#     check if pep includes flexible residues
    flex_pep_a,flex_pep_b,rigid_pep_a,rigid_pep_b = False,False,False,False
    nupeprigid=0
    if "pep" in segid_a:
        pep_unit = int(segid_a.lstrip("pep"))
        if int(resid_a) in flex[pep_unit-1]: #flexible residue
            flex_resid_a = flex[pep_unit-1].index(resid_a) +1
            flex_segid_a = "flx" + str(pep_unit)
            flex_pep_a = True
        else: # non-flex residue: retrieve coors
            rigid_pep_a=True
            nupeprigid+=1
            a_x, a_y, a_z = get_coor(atoms,pep_unit,resid_a,type_a)
            if a_x == None or a_y == None or a_z == None:
                return
    if "pep" in segid_b:
        pep_unit = int(segid_b.lstrip("pep"))
        if int(resid_b) in flex[pep_unit-1]: #flexible residue
            flex_resid_b = flex[pep_unit-1].index(resid_b) + 1
            flex_segid_b = "flx" + str(pep_unit)
            flex_pep_b = True
        else: # non-flex residue: retrieve coors
            rigid_pep_b=True
            nupeprigid+=1
            b_x, b_y, b_z = get_coor(atoms,pep_unit,resid_b,type_b)
            if b_x == None or b_y == None or b_z == None:
                return

    if nupeprigid==2: # noe only when pep
        line = "if @in1 .eq. pep then\n"
        line += simple_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k)
        line += "endif\n"    
    elif nupeprigid==1:
        line = "if @in1 .eq. pep then\n"
        line += simple_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k)
        line += "else\n"
        if flex_pep_a and rigid_pep_b:
            line += simple_mk_pnoe(b_x, b_y, b_z, flex_segid_a, flex_resid_a, type_a, rmin, rmax, k)
        elif flex_pep_b and rigid_pep_a:
            line += simple_mk_pnoe(a_x, a_y, a_z, flex_segid_b, flex_resid_b, type_b, rmin, rmax, k)
        elif (not flex_pep_a) and rigid_pep_b: # rigid prot to lig/cof
            line += simple_mk_pnoe(b_x, b_y, b_z, segid_a, resid_a, type_a, rmin, rmax, k)        
        elif rigid_pep_a and (not flex_pep_b): # rigid prot to lig/cof
            line += simple_mk_pnoe(a_x, a_y, a_z, segid_b, resid_b, type_b, rmin, rmax, k)        
        line += "endif\n"
        #write regular pnoe line when flex but when pep write noe 
    elif nupeprigid==0: #at least one flex pep
        line = "if @in1 .eq. pep then\n"
        line += simple_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k)
        line += "else\n"
        if flex_pep_a and (not flex_pep_b):
            line += simple_mk_noe(flex_segid_a, flex_resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k)
        elif flex_pep_b and (not flex_pep_a):
            line += simple_mk_noe(segid_a, resid_a, type_a, flex_segid_b, flex_resid_b, type_b, rmin, rmax, k)
        elif flex_pep_a and flex_pep_b:
            line += simple_mk_noe(flex_segid_a, flex_resid_a, type_a, flex_segid_b, flex_resid_b, type_b, rmin, rmax, k)
        line += "endif\n"
    return line

def mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k, flex, atoms):
    if ("pep" in segid_a.lower()) or ("pep" in segid_b.lower()):
        line = pep_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k, flex, atoms)
        if line == None:
            return
    else:
        line = simple_mk_noe(segid_a, resid_a, type_a, segid_b, resid_b, type_b, rmin, rmax, k)
    return line

def write_lines(df, fname):
    title = "\n".join(df["tot_line"])
    with open(fname,'a') as f:
        f.write(title)

def write_file_header(fname):
    with open(fname,'w') as f:
        title = """*************************************************************
*                Start userrestraints.str                   *
*************************************************************
* This file holds the noe restraints that the user supplied in
* the csv files named ./stream/consdef/userrestraints.csv
* it also contains some other mandatory restraints
* 

! NOTICE:
!  Restraints here are only read during docking.
!  They are NOT read during ligand preparation.
!

if @numcofactor .gt. 0 then
   cons harm force @kcf sele cofactors .and. .not. hydrogen show end
endif

if @water .eq. true then
   cons harm force @kwat sele segid CWAT .or. segid WTIN show end
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  NOE section

! do not reset since consensus is also noe
!noe
!   reset
!end

"""
        f.write(title)

def write_file_footer(fname):
    with open(fname,'a') as f:
        title = """!**********************************************************
!                  End userrestraints.str                 *
!**********************************************************

return

"""
        f.write(title)

def main():
    write_header()
    idir="../stream/consdef/"
    fname="userrestraints.csv"
    csvfname=idir + fname
    if path.isfile(csvfname):
        numligands, flexstr, proteinname = handle_input()
        flex, proteinunit = process_flexstr(flexstr)
        # process csv
        df, atoms = read_csv(csvfname,numligands,proteinname, proteinunit)
        # call_pep = df['col'].str.contains('pep').any() #bool - but read all protfiles?
        df["noe_header"] = df.apply(lambda x: "! noe " + str(x.name+2) + "\n", axis=1)
        df["ligsline1"] = df.apply(lambda x: mk_condition(x.ligands1), axis=1)
        df["line"] = df.apply(lambda x: mk_noe(x.segid_a, x.resid_a, x.type_a,
                                               x.segid_b, x.resid_b, x.type_b,
                                               x.rmin, x.rmax, x.k, flex, atoms), axis=1)
        if df['line'].isnull().values.any():
            quit()
        df["ligsline2"] = df.apply(lambda x: "endif\n" if len(x.ligands1) > 0 else "", axis=1)
        df["sele_warn"] = df.apply(lambda x: raise_sele_err(str(x.name+2)), axis=1)
        df["tot_line"] = df.apply(lambda x: x.noe_header + x.ligsline1 + x.line + x.ligsline2 + x.sele_warn, axis=1)
        # write file
        fname = "../stream/consdef/userrestraints.str"
        if path.isfile(fname):
            subprocess.call(["/usr/bin/chmod","+w",fname])
        write_file_header(fname)
        write_lines(df, fname)
        write_file_footer(fname)
    else:
        with open("enzydock.log",'a') as f:
            f.write("csv restraints file "+str(csvfname)+" does not exists. skipping...\n")

    fname=idir+"userrestraints.str"
    if path.isfile(fname):
        subprocess.call(["/usr/bin/chmod","-w",fname])
    write_footer()

if __name__ == "__main__":
    main()
 
