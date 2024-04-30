#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9

# Find MCS between ligands
# Can order the ligands 

#There are 4 cases:
#
#    1. known order - all the ligands compare to the first, which is the referance
#    2. known order - we compare every ligand with the previous ligand
#    3. no order - we find the referance, and compare all the ligands to the referance we have found (1)
#    4. no order - we find the best path, and compare every ligand with the previous ligand (2). 
#       here we can decide what the first ligand going to be and based on that we find the best path.

# If user provides csv files with matching pairs of atoms no MCS search os performed
#    file name: "../stream/consdef/user_consensus.csv"
#    Possible options from above: 1,2 (no reordering)
#    Structure of csv file: see in check_csv()

import os
import sys
import subprocess
import glob
import shutil
from rdkit import Chem
from rdkit.Chem import rdFMCS
from itertools import permutations
from read_write_pdb import PDB
import pandas as pd

def str2bool(var):
  var = var.lower()
  if var == "true":
    var = True
  else:
    var = False
  return var

def handle_input():
  # Handle input arguments and return values
  # set 0 @fpknownorder   <bool>
  # set 1 @fpmode         ("seq"/"seed")
  # set 2 @numligands     <int>
  # set 3 @currligand     <int>
  # set 4 @fromexp        <bool>
  # set 5 @fpseed         <int>; 0 for unknown seed
  # set 6 @anyele         <bool>
  # set 7 @runid          <str>
  # set 8 @iter           <int>
  # !set 1 if not yet defined
  # set 9 @multseed       <bool>
  input_arg = sys.argv
  known_order = input_arg[1].lower() # CHARMM sends in upper case parameters, change to lower
  mode = input_arg[2].lower() # CHARMM sends in upper case parameters, change to lower
  nlig = int(input_arg[3])
  currlig = int(input_arg[4])
  fromexp = input_arg[5].lower() # CHARMM sends in upper case parameters, change to lower
  fpseed = int(input_arg[6])
  anyele = input_arg[7].lower()
  runid = input_arg[8].lower()
  ligiter = int(input_arg[9])
  multseed = input_arg[10]

  known_order = str2bool(known_order)
  fromexp = str2bool(fromexp)
  anyele = str2bool(anyele)
  multseed = str2bool(multseed)

  return known_order, mode, nlig, currlig, fromexp, fpseed, anyele, runid, ligiter, multseed

def get_ref_coor(currlig, runid, ligiter):
  results_dir = "../results/"
  cluster_dir = results_dir + currlig + "_" + runid + "_" + "cluster" + str(ligiter) + "/"
  minfile = cluster_dir + "min_clust" + str(ligiter) + "_*.pdb"
  minfile = glob.glob(minfile)[0]
  return minfile

def write_log(line):
    print(line)
    # write to enzydock.log file
    with open("enzydock.log",'a') as f:
        f.write(line+"\n")

def check_csv():
  idir="../stream/consdef/"
  fname="user_consensus.csv"
  csvfname=idir + fname
  if os.path.isfile(csvfname):
    df = pd.read_csv(csvfname)
  else:
    df = None
  # structure of df:
  # ref1,lig1,ref2,lig2,...ref_n,lig_n
  # if no matching just leave empty cell
  # for each pair of columns, list pairs of atoms.
  # length of columns may vary between ligands but must be the same for ref_i, lig_i
  return df, csvfname

def create_mapping_dict(df, csvfname, nlig):
  ref_columns = [col for col in df.columns if col.startswith('ref')]
  lig_columns = [col for col in df.columns if col.startswith('lig')]
  if len(ref_columns) + len(lig_columns) != len(df.columns):
    line = f"redundant columns in {csvfname}: should only contain pairs of refn lign\n"
    line += "User's MCS restraints will not apply..."
    write_log(line)
    raise()
  if len(ref_columns) != len(lig_columns):
    line = f"missing columns in {csvfname}: should contain pairs of refb lign\n"
    line += "User's MCS restraints will not apply..."
    write_log(line)
    raise()
  # assuming correct structure if df and correct ordering:
  
  mapping_dicts = {}
  for ref_col, lig_col in zip(ref_columns, lig_columns):
    ref_values = df[ref_col].dropna().values
    lig_values = df[lig_col].dropna().values
    
    currref = int(ref_col.replace("ref",""))
    currligand = int(lig_col.replace("lig",""))
    if currref != currligand:
      line = f"missmatching columns in {csvfname}: {ref_col}, {lig_col}."
      line += "The CSV should contain pairs of refb lign"
      line += "User's MCS restraints will not apply..."
      write_log(line)
      raise()
    #idict = dict(zip(ref_values, lig_values))
    idict = {key.lower(): value.lower() for key, value in zip(ref_values, lig_values)}
    mapping_dicts[currligand] = idict
  for i in range(1, nlig+1):
    if not i in mapping_dicts:
      mapping_dicts[i]  ={}
  return mapping_dicts

def write_header():
    print("\nMCS ligand alignment program for EnzyDock\nCopyright © 2020 Dan T. Major\n")

def write_footer():
    print("\nEnd MCS ligand alignment program, returning to EnzyDock main\n")

def rename_ligand(ligs_dict):
  for org_l, new_l in ligs_dict.items():
  # ligs_dict = {1:1, 2:3, 3:2}
  
  #nligs=len(ligs_dict)
  #files=["../pdb/ligand_"+str(i)+".pdb" for i in range(1,nligs+1)]
  #files+=["../local_top/ligand_"+str(i)+".str" for i in range(1,nligs+1)]

    """    
    @pdbDIR/ligand_@currligand.pdb #ligand
    @topDIR/ligand_@currligand.str #ligand
    
    @topDIR/ligand_@currligand_@cofact@@d_patch.str #ligand-cof
    @stpDIR/cofact@d_patch_@currligand.str          #ligand-cof
    
    @stpDIR/cofact@c_patch_ng_@currligand.str #cof-prot
    @stpDIR/cofact@c_patch_g_@currligand.str  #cof-prot
                                              # + definitions somewhere in local_top - not relevant
    
    @stpDIR/patch_ng_@currligand.str #prot
    @stpDIR/patch_g_@currligand.str  #prot
                                     # + definitions in set_prot_patch.str !!
                                     # + definitions in local_top - not relevant
    
    @stdDIR/@currligand_consensus_data.str
    @stdDIR/@currligand_fp_data.str
    
    ?? What to do with covalent docking #lig-prot --> add note in param_set_check (dictate it as the first)
    ?? free-style userrestraints - block option to MCS re-ordering --> add warning to user
    
    update following files:
    |- userestraints.csv
    |- userparam.str
       |- set lastinter @numligands
       |- set qmcharge1 0
       |- set qmcharge2 0    ! ligand
       |- set qmcharge3 0    ! ligand
    |- param_set_check.str
    |- set_prot_patch.str
    """
    
  # 2 loops: move all to original; copy+rename - from original to avoid overwriting
    #TODO: add ligand rtf/prm files
    shutil.move("../pdb/ligand_"+str(org_l)+".pdb", "../pdb/ligand_"+str(org_l)+"_org.pdb")
    shutil.move("../local_top/ligand_"+str(org_l)+".str", "../local_top/ligand_"+str(org_l)+"_org.str")
    
    #f="../stream/consdef/"+str(org_l)+"_consensus_data.str"
    #f="../stream/consdef/"+str(org_l)+"_fp_data.str"
    
    fs=glob.glob("../stream/patching/patch_*g_"+str(org_l)+".str")
    fs+=glob.glob("../local_top/ligand_"+str(org_l)+"_*_patch.str")
    fs+=glob.glob("../stream/patching/cofact*_patch_ng_"+str(org_l)+".str")
    fs+=glob.glob("../stream/patching/cofact*_patch_g_"+str(org_l)+".str")
    fs+=glob.glob("../stream/patching/cofact*_patch_"+str(org_l)+".str")
    for f0 in fs:
      f1=f0.replace(".str","_org.str")
      #if os.path.exists(f0):
      shutil.move(f0,f1)

  for org_l, new_l in ligs_dict.items():
    #f="../pdb/ligand_"+str(org_l)+".pdb"
    shutil.copy2("../pdb/ligand_"+str(org_l)+"_org.pdb", "../pdb/ligand_"+str(new_l)+".pdb")
    shutil.copy2("../local_top/ligand_"+str(org_l)+"_org.str", "../local_top/ligand_"+str(new_l)+".str")
    
    #f="../stream/consdef/"+str(org_l)+"_consensus_data.str"
    #f="../stream/consdef/"+str(org_l)+"_fp_data.str"
      
    fs=glob.glob("../local_top/ligand_"+str(org_l)+"_*_patch_org.str")
    for f0 in fs:
      f1=f0.replace("../local_top/ligand_"+str(org_l)+"_","../local_top/ligand_"+str(new_l)+"_")
      f1=f1.replace("_patch_org.str","_patch.str")
      shutil.copy2(f0,f1)

    fs=glob.glob("../stream/patching/*patch_*"+str(org_l)+"_org.str")
    for f0 in fs:
      f1=f0.replace(str(org_l)+"_org.str",str(new_l)+".str")
      shutil.copy2(f0,f1)

    # fs=glob.glob("../stream/patching/patch_*g_"+str(org_l)+"_org.str")
    # for f0 in fs:
      # f1=f0.replace("g_"+str(org_l)+"_org.str","g_"+str(new_l)+".str")
      # shutil.copy2(f0,f1)

    # fs=glob.glob("../stream/patching/cofact*_patch_*"+str(org_l)+"_org.str")
    # for f0 in fs:
      # f1=f0.replace(str(org_l)+"_org.str",str(new_l)+".str")
      # shutil.copy2(f0,f1)

    # fs=glob.glob("../stream/patching/cofact*_patch_g_"+str(org_l)+".str")
    # for f0 in fs:
      # f1=f0.replace("_patch_g_"+str(org_l)+".str","_patch_g_"+str(new_l)+".str")
      # shutil.copy2(f0,f1)
    # fs=glob.glob("../stream/patching/cofact*_patch_"+str(org_l)+".str")
    # for f0 in fs:
      # f1=f0.replace("_patch_"+str(org_l)+".str","_patch_"+str(new_l)+".str")
      # shutil.copy2(f0,f1)
        
def handle_noes(ligs_dict):
  fname0="../stream/consdef/userrestraints.csv"
  if os.path.exists(fname0):
    fname1="../stream/consdef/userrestraints.str"
    with open(fname1) as f:
      lines=f.readlines()
    newlines=""
    for l in lines:
      if l.startswith("if @currligand .eq. "):
        old_l=l.lstrip("if @currligand .eq. ").rstrip(" set donoe true\n")
        new_l=ligs_dict[int(old_l)]
        l1="if @currligand .eq. "+str(new_l)+" set donoe true\n"
        newlines+=l1
      else:
        newlines+=l#+"\n"
    
    shutil.move(fname1,fname1.replace(".str","_org.str"))
    
    with open(fname1,"w") as f:
      f.write(newlines)

# def handle_ppat(ligs_dict):
  # fname0="../stream/patching/set_prot_patch.str"
  # with open(fname0) as f:
    # lines=f.readlines()
  # newlines=""
  # check=True
  # for l in lines:
    # if check:
      # if "loop over ligands" in l:
        # check=False
        # newlines+=l#+"\n"
      # else:
        # if l.isspace() or (l.startswith("!")) or (l.startswith("*")):
          # newlines+=l#+"\n"
        # elif "set patchcountligp" in l:
          # old_l=l.split()[1].lstrip("patchcountligp")[:-1]
          # new_l=ligs_dict[int(old_l)]
          # l1=l.replace("patchcountligp"+str(old_l),"patchcountligp"+str(new_l))
          # newlines+=l1#+"\n"
        # elif "set patchnamelig" in l:
          # old_l=l.split()[1].lstrip("patchnamelig").split("p")[0]
          # new_l=ligs_dict[int(old_l)]
          # l1=l.replace("patchnamelig"+str(old_l),"patchnamelig"+str(new_l))
          # newlines+=l1#+"\n"
        # elif "set patchreslig" in l:
          # old_l=l.split()[1].lstrip("patchreslig").split("p")[0]
          # new_l=ligs_dict[int(old_l)]
          # l1=l.replace("patchreslig"+str(old_l),"patchreslig"+str(new_l))
          # newlines+=l1#+"\n"
    # else:
      # newlines+=l#+"\n"
  
  # shutil.move(fname0,fname0.replace(".str","_org.str"))
  
  # with open(fname0,"w") as f:
    # f.write(newlines)

def handle_ppat(ligs_dict):
  fname0="../stream/patching/set_prot_patch.str"
  prelines="""* Temporary file to define patches variables
* To be merged into userparam.str and param_set_check.str
* RS 13/05/2021
* 

! give default values for the count - make all of them 0
! then re-define them

! loop over ligands
set l 1 ! lign 1
label pre_countp_lig
    ! loop over protein units for a specific ligand
    set prot 1 ! p 1
    label pre_countp_prot
        set patchcountligp@@{l}@@{prot} 0
        incr prot
    if @prot .le. @proteinunit goto pre_countp_prot
    incr l
if @l .le. @numligands goto pre_countp_lig

"""
  with open(fname0) as f:
    lines=f.readlines()
  newlines=prelines
  check=True
  for l in lines:
    if "loop over ligands" in l:
      break
    if l.isspace() or (l.startswith("!")) or (l.startswith("*")):
      continue
    elif "set patchcountligp" in l:
      old_l=l.split()[1].lstrip("patchcountligp")[:-1]
      new_l=ligs_dict[int(old_l)]
      l1=l.replace("patchcountligp"+str(old_l),"patchcountligp"+str(new_l))
      newlines+=l1#+"\n"
    elif "set patchnamelig" in l:
      old_l=l.split()[1].lstrip("patchnamelig").split("p")[0]
      new_l=ligs_dict[int(old_l)]
      l1=l.replace("patchnamelig"+str(old_l),"patchnamelig"+str(new_l))
      newlines+=l1#+"\n"
    elif "set patchreslig" in l:
      old_l=l.split()[1].lstrip("patchreslig").split("p")[0]
      new_l=ligs_dict[int(old_l)]
      l1=l.replace("patchreslig"+str(old_l),"patchreslig"+str(new_l))
      newlines+=l1#+"\n"
  
  endlines="""!!!!!!!!!!!!!!!!!!!
! param_set_check.str

! give default values for the count
! if patch is defined but not count, it'll be ignored
! additional checks might be needed

! proteinunit is known
! numligands is known

! loop over ligands
set l 1 ! lign 1
label countp_lig
    ! loop over protein units for a specific ligand
    set prot 1 ! p 1
    label countp_prot
        set test @patchcountligp@@{l}@@{prot}
        if @test .gt. 0 then
            ! set felx residue (different indexing)
            set respatch 1
            ! loop over patches for a specific unit and a specific ligand
            label patch_set_flex
                goto assign_flex
                label after_assign_flex
                !!!! DEBUG PRINT STARTS
                echo @patchflexlig@@{l}p@@{prot}n@@{respatch}
                echo @patchnamelig@@{l}p@@{prot}n@@{respatch}
                !!!! DEBUG PRINT ENDS
                incr respatch
            if respatch .le. @patchcountligp@@{l}@@{prot} goto patch_set_flex
        endif
        incr prot
    if @prot .le. @proteinunit goto countp_prot
    incr l
if @l .le. @numligands goto countp_lig

return
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the index of the residue among the flexible residues of the protein unit
!!! patch definitions (from above)
! ligand : @l
! protein unit : @prot
! patch index [serial number] : @respatch
! residue number: @{patchreslig@@{l}p@@{prot}n@@{respatch}}

!!! flexible residues definitions (from above)
! number of flexible residues in relevant protein unit: @numflex@@prot
! residue number of the Nth flexible residue : @flex@@prot@@n

label assign_flex
! make sure patch name and residue exists
set test 0
set test @?patchnamelig@@{l}p@@{prot}n@@{respatch}
if @test .eq. 0 then
    echo EnzyDock WARNING: No patch name found for ligand @l    
    echo EnzyDock WARNING: protein unit @prot patch number @respatch
    ! Update readable log file
!    open append unit 33 form name enzydock.log
    write title unit 33
* EnzyDock WARNING: No patch name found for ligand @l
* EnzyDock WARNING: protein unit @prot patch number @respatch
* Terminating EnzyDock run...
*
    close unit 33
    stop
endif

set test 0
set test @?patchreslig@@{l}p@@{prot}n@@{respatch}
if @test .eq. 0 then
    echo EnzyDock WARNING: No residue number found for ligand @l    
    echo EnzyDock WARNING: protein unit @prot patch number @respatch
    ! Update readable log file
!    open append unit 33 form name enzydock.log
    write title unit 33
* EnzyDock WARNING: No residue number found for ligand @l
* EnzyDock WARNING: protein unit @prot patch number @respatch
* Terminating EnzyDock run...
*
    close unit 33
    stop
endif

set assigned false
set flexind 1
! loop over flexible residues of a specific ligand and protein unit
label match_flex
    if @patchreslig@@{l}p@@{prot}n@@{respatch} .eq. @flex@@{prot}@@{flexind} then
        set patchflexlig@@{l}p@@{prot}n@@{respatch} @flexind
        set patchflexnamel@@{l}p@@{prot}i@@{flexind} @patchnamelig@@{l}p@@{prot}n@@{respatch}
        set assigned true
    endif
    incr flexind
if @flexind .le. @numflex@@{prot} goto match_flex

! make sure flexible residue was found

if @assigned .eq. true then
    goto after_assign_flex
else
    echo EnzyDock WARNING: Residue that needs a patch is not defined as flexible:
    echo EnzyDock WARNING: Ligand state: @l; Protein unit: @prot; Patch nuber: @respatch
    echo EnzyDock WARNING: Residue number: @patchreslig@@{l}p@@{prot}n@@{respatch}
    ! Update readable log file
!    open append unit 33 form name enzydock.log
    write title unit 33
* EnzyDock WARNING: Residue that needs a patch is not defined as flexible:
* EnzyDock WARNING: Ligand state: @l; Protein unit: @prot; Patch nuber: @respatch
* EnzyDock WARNING: Residue number: @patchreslig@@{l}p@@{prot}n@@{respatch}
* Terminating EnzyDock run...
*
    close unit 33
    stop
endif

"""
  newlines+=endlines
  shutil.move(fname0,fname0.replace(".str","_org.str"))
  
  with open(fname0,"w") as f:
    f.write(newlines)

def qm_charges(ligs_dict):
  fname0="../stream/userparam.str"
  with open(fname0) as f:
    lines=f.readlines()
  newlines=""
  for l in lines:
    if ("qmcharge" in l) and ("set" in l):
      if l.find("!") > -1:
        if l.find ("!") < l.find("set"):
          continue
      wrds=l.split()
      old_l=wrds[wrds.index("set")+1].lstrip("qmcharge")
      new_l=ligs_dict[int(old_l)]
      l1=l.replace("qmcharge"+str(old_l),"qmcharge"+str(new_l))
      newlines+=l1#+"\n"
  fname1="../stream/patching/set_prot_patch.str"
  with open(fname1,"a+") as f:
    f.write(newlines)

def reorder_ligs(ligs_dict):
  with open("enzydock.log",'a') as f:
     f.write("MCS program reorders ligands:...\n")
     for old,new in ligs_dict.items():
       f.write("ligand_"+str(old)+" --> ligand_"+str(new)+"\n")
  rename_ligand(ligs_dict)
  handle_noes(ligs_dict)
  handle_ppat(ligs_dict)
  qm_charges(ligs_dict)

def charmm2rdkit(pdbf):
    # First identify unused backup filename 
    tmpfilename = 'tmp'
    if (os.path.exists(tmpfilename)):
       tmpfilename0 = tmpfilename
       tmpfilename = tmpfilename0 + '_0'
       if (not os.path.exists(tmpfilename)):
          exist = False
       else:
          exist = True
          i = 1
       while (exist):
          tmpfilename = tmpfilename0 + '_' + str(i)
          if (not os.path.exists(tmpfilename)):
             exist = False
          else:
             i += 1

    pdb_file = pdbf
    subprocess.call(["/bin/rm","-f",tmpfilename])
    pdb = PDB(pdb_file, tmpfilename)
    atoms = pdb.get_atoms(to_dict=False)
    pdb.add_element()
    pdb.write_pdb(tmpfilename)
    subprocess.call(["/bin/mv",tmpfilename,pdb_file])

#BASIC FUNCTIONS FOR THE CASES 

#Calculate the common structure of the two molecules (k)
def mcs_calculate(a, b, anyele):
    k = rdFMCS.FindMCS([a, b], bondCompare=rdFMCS.BondCompare.CompareAny)
    if anyele:
        k = rdFMCS.FindMCS([a, b], bondCompare=rdFMCS.BondCompare.CompareAny,
                           atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom)
    return k
    
def atom_name(a, b, k):
    #Get tuple of indices of the atoms from the REF MOLECULE & TEST MOLECULE
    a_indices = (a.GetSubstructMatches(Chem.MolFromSmarts(k), uniquify = False, useChirality=True))[0]
    b_indices = (b.GetSubstructMatches(Chem.MolFromSmarts(k), uniquify = False, useChirality=True))[0]

    #For a given index, we get the atom's name according to PDB, and it's coordinates
    a_atom_names = []
    b_atom_names = []

    for i in a_indices:
        a_atom_names.append(a.GetAtomWithIdx(i).GetMonomerInfo().GetName())
    for i in b_indices:
        b_atom_names.append(b.GetAtomWithIdx(i).GetMonomerInfo().GetName())
        
    return a_atom_names, b_atom_names, a_indices, b_indices

#Get the coordinates of the atoms in a given ligand
def atom_coordi(a_indices, b_indices, a_coor):
    a_atoms_xyz = []
    b_atoms_xyz = []    
    for i in a_indices:
        positions = a_coor.GetConformer().GetAtomPosition(i)
        a_atoms_xyz.append([positions.x, positions.y, positions.z])

    return a_atoms_xyz

def number_of_atoms(a_atoms_xyz):
    return len(a_atoms_xyz)

#Get a file of the characters of the reference and the ligand
def writefunc(a_atom_names, a_atoms_xyz, b_atom_names, num_lig):
    fname = "../stream/consdef/"+str(num_lig) + "_fp_data.str"
    line = "* Fingerprint data for docking definitions" + "\n"
    line += "* Will be written by python code" + "\n"
    line += "* fpatom is finperprint atom" + "\n"
    line += "* latom is ligand atom" + "\n"
    #line += "* Total i in number of latom@i@j must match number of ligands" + "\n"
    line += "* \n" + "\n"
    line += "set numrefatom " + str(len(a_atom_names)) + "\n"
    a = 1
    for i in range(len(a_atom_names)):
        line += "set fpatom" + str(a) + " " + str(a_atom_names[i])+ "\n"
        line += "set xfpatom" + str(a) + " "+ str(a_atoms_xyz[i][0])+ "\n"
        line += "set yfpatom" + str(a) + " " + str(a_atoms_xyz[i][1])+ "\n"
        line += "set zfpatom" + str(a) + " " + str(a_atoms_xyz[i][2])+ "\n"
        line += "set latom" + str(a) + " " + str(b_atom_names[i])+ "\n" +"\n"
        a += 1
    with open (fname, "w") as f:
        f.write(line)
        
def calc_properties(a, a_coor, b, num_lig, anyele):
    #a, b, a_coor = get_files(reference_struct, reference_coor, ligand)
    k = mcs_calculate(a, b, anyele)
    k = k.smartsString
    a_atom_names, b_atom_names, a_indices, b_indices = atom_name(a, b, k)
    a_atoms_xyz = atom_coordi(a_indices, b_indices, a_coor)
    writefunc(a_atom_names, a_atoms_xyz, b_atom_names, num_lig)

def names_to_xyz(pdbf):
    with open(pdbf) as f:
      lines = f.readlines()
    pdb = {}
    for line in lines:
      if line.startswith('ATOM') or line.startswith('HETATM'):
        #atname = line[11:16].strip()
        atname = line[11:16].strip().lower()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        pdb[atname] = {"x":x, "y":y, "z":z}
    return pdb

def dict2lists(currdict, ref_xyz):
    a_atom_names = []
    a_atoms_xyz = []
    b_atom_names = []
    for ref,lig in currdict.items():
        if ref in ref_xyz:
            coor = ref_xyz[ref]
            a_atom_names.append(ref)
            a_atoms_xyz.append([coor["x"], coor["y"], coor["z"]])
            b_atom_names.append(lig)
        else:
            line = f"atom {ref} was requested but not found in reference pdb file"
            write_log(line)
    return a_atom_names, a_atoms_xyz, b_atom_names

def u_calc_properties(currligand, currdict, reference_coor):
    ref_xyz = names_to_xyz(reference_coor)
    a_atom_names, a_atoms_xyz, b_atom_names = dict2lists(currdict, ref_xyz)
    writefunc(a_atom_names, a_atoms_xyz, b_atom_names, currligand)
    
def MCS_matrix(number_of_ligands, from_exp, anyele):
#This function has a call to files - pay attention when switching to Enzydock
    l = number_of_ligands
    ligs = []
    common_atoms_matrix = [] 
    #make a list from the ligands that we insert to the code
    if from_exp:
        i=0
        ligand = "../pdb/ligand_" + str(i) + ".pdb"
        ligs.append(ligand)        
    for i in range(1, int(l)+1):
        ligand = "../pdb/ligand_" + str(i) + ".pdb"
        ligs.append(ligand)
    for i in range(len(ligs)):
        mol_object = Chem.MolFromPDBFile(ligs[i])
        ligs[i] = mol_object

    #find the reference - the ligand that have the max common atoms & bonds with the others
    for a, i in enumerate(ligs):
        ks = []
        for b, j in enumerate(ligs):
            if a < b:
                k = mcs_calculate(i, j, anyele)
                k = k.numAtoms
            elif a == b:
                k = 0
            else:
                k = common_atoms_matrix[b][a]
            ks.append(k)
        common_atoms_matrix.append(ks)

    return common_atoms_matrix

### case 1: known order - all the ligands compare to the first, which is the reference
def all_with_one(l, from_exp, anyele, runid, ligiter, multseed):
#This function has a call to files - pay attention when switching to Enzydock
    if from_exp == True:  #(the reference is 0)
        reference_struct = Chem.MolFromPDBFile("../pdb/ligand_0.pdb")
        # same if multseed==True
        reference_coor = Chem.MolFromPDBFile("../pdb/ligand_0.pdb")
        print("MCS reference: ../pdb/ligand_0.pdb")
        for i in range(1, int(l)+1):
            i = str(i)
            num_lig = i
            ligand = Chem.MolFromPDBFile("../pdb/ligand_" + i + ".pdb")
            calc_properties(reference_struct, reference_coor, ligand, num_lig, anyele)
    elif from_exp == False: #(the reference is 1)
        reference_struct = Chem.MolFromPDBFile("../pdb/ligand_1.pdb")
        ref_pdb="../results/lig_mindock_1.pdb"
        if multseed:
            ref_pdb=get_ref_coor("1", runid, ligiter)
        print(f"MCS reference: {ref_pdb}")
        charmm2rdkit(ref_pdb)
        reference_coor = Chem.MolFromPDBFile(ref_pdb)
        for i in range(2, l+1):
            i = str(i)
            num_lig = i
            ligand = Chem.MolFromPDBFile("../pdb/ligand_" + i + ".pdb")        
            calc_properties(reference_struct, reference_coor, ligand, num_lig, anyele)

def u_all_with_one(nlig, from_exp, mapping_dicts, runid, ligiter, multseed):
    if from_exp == True:  #(the reference is 0)
        # same if multseed==True
        reference_coor = "../pdb/ligand_0.pdb"
        start = 1
    elif from_exp == False: #(the reference is 1)
        reference_coor = "../results/lig_mindock_1.pdb"
        if multseed:
            reference_coor = get_ref_coor("1", runid, ligiter)
        start = 2
    print(f"MCS reference (user's matching): {reference_coor}")
    for i in range(start, int(nlig)+1):
        if not (i in mapping_dicts):
            continue
        currdict = mapping_dicts[i]
        u_calc_properties(i, currdict, reference_coor)

### case 2: known order - we compare every ligand with the previous ligand
def ligand_with_previous(num_lig, anyele, runid, ligiter, multseed):
    if num_lig == 1:
        # same if multseed==True
        print("MCS reference: ../pdb/ligand_0.pdb")
        reference_coor = Chem.MolFromPDBFile("../pdb/ligand_0.pdb")
    else:
        ref_pdb="../results/lig_mindock_" + str(num_lig-1) + ".pdb"
        if multseed:
            ref_pdb=get_ref_coor(str(num_lig-1), runid, ligiter)
        print(f"MCS reference: {ref_pdb}")
        charmm2rdkit(ref_pdb)
        reference_coor = Chem.MolFromPDBFile(ref_pdb)
    
    ligand = Chem.MolFromPDBFile("../pdb/ligand_" + str(num_lig) + ".pdb")
    reference_struct = Chem.MolFromPDBFile("../pdb/ligand_" + str(num_lig-1) + ".pdb")
    calc_properties(reference_struct, reference_coor, ligand, num_lig, anyele)

def u_ligand_with_previous(num_lig, mapping_dicts, runid, ligiter, multseed):
    if num_lig == 1:
        # same if multseed==True
        reference_coor = "../pdb/ligand_0.pdb"
    else:
        reference_coor = "../results/lig_mindock_" + str(num_lig-1) + ".pdb"
        if multseed:
            reference_coor = get_ref_coor(str(num_lig-1), runid, ligiter)
    if not (num_lig in mapping_dicts):
        line=f"No restraints found for ligand {num_lig}, skipping restraints..."
        write_log(line)
    else:
        currdict = mapping_dicts[num_lig]
        print(f"MCS reference (user's matching): {reference_coor}")
        u_calc_properties(num_lig, currdict, reference_coor)

### case 3: no order - we find the reference, and compare all the ligands to the reference we have found  
def find_ref_in_ligands(number_of_ligands, anyele):
#This function has a call to files - pay attention when switching to Enzydock
    common_atoms_matrix = MCS_matrix(number_of_ligands, False, anyele)
    ref = max(range(len(common_atoms_matrix)), key=lambda i: sum(common_atoms_matrix[i]))
    ref += 1
    
    #for i in range(1,1+number_of_ligands):
    
    ##now that we found the ref, we need to call it "ligand num 1"
    #if ref != 1:
    #    os.rename('ligand_1.pdb', 'ligand_temporary.pdb')
    #    os.rename('ligand_' + str(ref) + '.pdb', 'ligand_1.pdb')
    #    os.rename('ligand_temporary.pdb', 'ligand_' + str(ref) + '.pdb')     
    return ref

### case 4: no order - we find the best path, and compare every ligand with the previous ligand.
def best_path(common_atoms_matrix, from_exp, chosen_head = 0):
    
    #IN THIS FUNCTION WE CREATE A LIST OF ALL THE COMBNATIONS THAT CAN BE CREATED FROM THE LIGANDS.
    #THEN, WE RANK EVERY COMBINATION WITH THE SUM OF THE K BETWEEN THE LIGANDS.
    #THE BEST RANK = THE BEST COMBINATION = THE BEST PATH.
    
    #chosen_head = the ligand that I CHOSE to be the first in the path. must be int!!!
    
    mat = []
    for i in range(len(common_atoms_matrix)):
        mat.append(i)

    perms = permutations(mat)

    ligands_nums = []
    for j in list(perms):
        ligands_nums.append(j)

    #we cut the list of the combinations by 2 (because there is  (0 1 2 3) and (3 2 1 0) [== same])

    ligands_nums_short = []

    if from_exp:
        new_list = [tup for tup in ligands_nums if tup[0] == 0]
    else:
        if chosen_head == 0: # no desired seed
            for i in ligands_nums:
                if i[-1] < i[0]:
                    ligands_nums_short.append(i)
            new_list = ligands_nums_short
        else:
            new_list = [tup for tup in ligands_nums if tup[0] == chosen_head-1] # if not from exp ligand 1 is in row #0


    ranks = []    
    for i in new_list: # i is a sequence of ligands
        counter = 0
        for a, b in zip(i[:-1], i[1:]):
            counter += common_atoms_matrix[a][b]
        ranks.append(counter)

    max_sum = ranks.index(max(ranks)) #find the maximal sum (one of its occurrences)

    best_ligands_path = new_list[max_sum] #goes to the list of the ligands combinations and takes the best

    #Now we want to check which path is better - the path we calculated or the inverse of it.
    #So if the numbers of common atoms are, for example (4, 5, 7, 10), 
    #we will want the inverse path because we want to start from the highest number of common atoms.

    if from_exp:
        best_ligands_path = best_ligands_path[1:]
    else: #add a "shift"
        if chosen_head ==  0: #no chosen head when not from exp, so select the direction 
            if common_atoms_matrix[best_ligands_path[0]][best_ligands_path[1]] < common_atoms_matrix[best_ligands_path[-1]][best_ligands_path[-2]]:
                best_ligands_path = best_ligands_path[::-1]
        
        # add shift
        best_ligands_path = tuple(i + 1 for i in best_ligands_path)

    return best_ligands_path

#def rename_ligands_files(best_ligands_path):
#    for indice, ligand in enumerate(best_ligands_path, 1):
#        #construct the old and new file names
#        old_file_name = 'ligand_' + str(ligand+1) + '.pdb'
#        new_file_name = 'ligand_' + str(ligand+1) + '_org.pdb'
#        #make a copy of the original file
#        shutil.copy(old_file_name, new_file_name)
#
#        os.rename('ligand_' + str(ligand+1) + '.pdb', 'ligand_' + str(indice) + '_new.pdb')
#    for indice in range(1, len(best_ligands_path)+1):
#        os.rename('ligand_' + str(indice) + '_new.pdb', 'ligand_' + str(indice) + '.pdb')

def make_ligs_dict(best_ligands_path):
    ligs_dict = {}
    for new_n, curr_n in enumerate(best_ligands_path, 1):
        if curr_n != 0:
            ligs_dict[curr_n] = new_n
    return ligs_dict

def ref_to_ligs_dict(ref,nlig):
    l = [ref]
    for i in range (2,nlig+1):
        if i != ref:
            l.append(i)
        else:
            l.append(1)
    return make_ligs_dict(l)

## TODO: in enzydock:
# re-read ../stream/patching/set_prot_patch.str !!
    
def main():
    write_header()
    known_order, mode, nlig, currlig, fromexp, fpseed, anyele, runid, ligiter, multseed = handle_input()
    user=False
    # check if user provided csv with pairs of atoms: (no reordering allowed)
    df, csvfname = check_csv()
    if not(df is None):
        user=True
        mapping_dicts = create_mapping_dict(df, csvfname, nlig)

    # selection of the case (known/unknown order; seq/seed(=all to 1))
    if known_order == True or user:
        if mode == "seed":
            if user:
                u_all_with_one(nlig, fromexp, mapping_dicts, runid, ligiter, multseed) #number of ligands
            else:
                all_with_one(nlig, fromexp, anyele, runid, ligiter, multseed) #number of ligands
        elif mode == "seq":
            if user:
                u_ligand_with_previous(currlig, mapping_dicts, runid, ligiter, multseed) #current ligand
            else:
                ligand_with_previous(currlig, anyele, runid, ligiter, multseed) #current ligand
    elif known_order == False:
        ## TODO: make sure enzydock calls the script again to actually create the reference files
        if mode == "seed":
            ref = find_ref_in_ligands(nlig, anyele) #switch between the ref and ligand files names
            ligs_dict = ref_to_ligs_dict(ref,nlig)
            reorder_ligs(ligs_dict)
        elif mode == "seq":
            common_atoms_matrix = MCS_matrix(nlig, fromexp, anyele)
            best_ligands_path = best_path(common_atoms_matrix, fromexp, chosen_head = fpseed)
            ligs_dict = make_ligs_dict(best_ligands_path)
            reorder_ligs(ligs_dict)
    write_footer()

if __name__ == "__main__":
    main()
