#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
# To test as standalone use:
# ./prepare.py /private/chem/silcsbio.2022.1/cgenff/cgenff
# 
# Copyright © 2022 Dan T. Major

import sys
import os
sys.path.append(".")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import pandas as pd
from openbabel import pybel

class Prep:

    def __init__(self, mol, name, pdb_dir, local_top):
        self.mol = mol
        self.name = str(name)
        self.pdb_dir = pdb_dir
        self.local_top = local_top
        self.mol_file = pdb_dir + name + '.mol'
        self.mol2_file = pdb_dir + name + '.mol2'
        self.str = local_top + name + '.str'
        self.pdb_file = pdb_dir + name + '.pdb'
        input_arg = sys.argv
        self.cgenff = input_arg[1].lower()

    def write_to_mol(self):
        AllChem.EmbedMolecule(self.mol)
        print(Chem.MolToMolBlock(self.mol), file=open(self.pdb_dir + self.name + '.mol', 'w+'))
        
    def label_atoms(self):
        for i, atom in enumerate(self.mol.GetAtoms()):
            atom.SetAtomMapNum(i + 1)
        return self

    def mol2_to_pdb(self):
        for mymol in pybel.readfile('mol2', self.pdb_dir + self.name  + '.mol2'):
            output = pybel.Outputfile('pdb', self.pdb_dir + self.name + '.pdb', overwrite=True)
            print(output.write(mymol))
            output.close()

    def mol_to_mol2(self):
        try:
            for mymol in pybel.readfile('mol', self.pdb_dir + self.name + '.mol'):
                output = pybel.Outputfile('mol2', self.pdb_dir + self.name + '.mol2', overwrite=True)
                print(output.write(mymol))
                output.close()
        except OSError:
            pass

    def run_cgenff(self):
        #cgenff = '/private/chem/silcsbio.2022.1/cgenff/cgenff ' # define in scripts/define_system_var.sh
        os.system(self.cgenff + ' ' + self.pdb_dir + self.name + '.mol2 1>' + self.local_top + self.name + '.str 2>>' +self.local_top + self.name +'.log')

    def rename_lig(self):
        with open(self.str , "r") as f:
            lines = f.readlines()
            f.close()
        index_line = [j for j in range(len(lines)) if 'RESI' in lines[j]][0]
        new_line = [lines[j].replace(lines[j].split()[1], 'LIG') for j in range(len(lines)) if 'RESI' in lines[j]][0]
        lines[index_line] = new_line
        f = open(self.str , "w")
        lines = "".join(lines)
        f.write(lines)
        f.close()

    def fix_mol2(self):
        with open(self.mol_file, "r") as f:
            lines = f.readlines()
            f.close()
        atom_numbering = [lines[j].split()[-3] for j in range(len(lines)) if (lines[j].split() and len(lines[j].split()) == 16)]

        with open(self.mol2_file, "r") as f:
            lines = f.readlines()
            f.close()
        new_lines = [lines[j].replace(lines[j].split()[1], lines[j].split()[1] + str(atom_numbering[j-7]),1)  if  (lines[j].split() and len(lines[j].split()) == 9) else lines[j] for j in range(len(lines))]
        f = open(self.mol2_file, "w")
        lines = "".join(new_lines)
        f.write(lines)
        f.close()

    def fix_pdb(self):
        with open(self.pdb_file, "r") as file:
            filedata = file.read()
        filedata = filedata.replace('UNL', 'LIG')
        with open(self.pdb_file, "w") as file:
            file.write(filedata)

    def delete_mol(self):
        os.remove(self.mol2_file)
        os.remove(self.mol_file)

def write_header():
    print("\nLigand prep program for EnzyDock\nCopyright © 2022 Dan T. Major\n")

def write_footer():
    print("\nEnd ligand prep program, returning to EnzyDock main\n")

def main():
    write_header()
    pdb_dir = '../pdb/'
    local_top = '../local_top/'
    stream = '../stream/'
    input_arg = sys.argv
    cgenff = input_arg[1].lower() 
    #f = open(stream + 'prepare.str', 'w')
    first_lines = '************************************************************  \n*                   prepare.str          *   \n************************************************************ \n*stop enzydock if prepare.py didnt work correctly \n*\n '

    try: 
        ligands = pd.read_csv(pdb_dir+'ligands.csv')
    except FileNotFoundError: 
        # RS: Add report to enzydock.log file
        line = "EnzyDock WARNING: ligands.csv file not found\n"
        line += "Terminating EnzyDock run...\n"
        with open("enzydock.log",'a') as f:
            f.write(line)
        # RS: End of report to enzydock.log file
        
        f = open(stream + 'prepare.str', 'w')
        f.write(first_lines +'echo csv file not found \n stop \n')
        f.close()
        #print('FileNotFoundError')
        #sys.exit()
    try:  
        ligands['mol'] = ligands['smile'].apply(lambda x: Chem.AddHs(Chem.MolFromSmiles(x)))
    except KeyError:
        # RS: Add report to enzydock.log file
        line = "EnzyDock WARNING: ligands.csv file in a wrong format:\n"
        line += "ligands.csv should include two columns: \"smile\" and \"name\"\n"
        line += "Terminating EnzyDock run...\n"
        with open("enzydock.log",'a') as f:
            f.write(line)
        # RS: End of report to enzydock.log file
        
        f = open(stream + 'prepare.str', 'w')
        f.write(first_lines+'echo the format of the csv file is not correct- should include one columns with name "smile" and second column "name"   \n stop \n')
        f.close()
        #print('keyerror') 
        #sys.exit() 
    ligands['m'] = ligands.apply(lambda x: Prep(x['mol'], 'ligand_' + str(x.name+1), pdb_dir, local_top), axis =1)
    ligands['m_l'] = ligands.apply(lambda x: x.m.label_atoms(), axis=1)
    ligands.apply(lambda x: x.m_l.write_to_mol(), axis=1)
    ligands.apply(lambda x: x.m.mol_to_mol2(), axis=1)
    ligands.apply(lambda x: x.m.fix_mol2(), axis=1)
    ligands.apply(lambda x: x.m.mol2_to_pdb(), axis=1)
    ligands.apply(lambda x: x.m.run_cgenff(), axis=1)
    ligands.apply(lambda x: x.m.rename_lig(), axis=1)
    ligands.apply(lambda x: x.m.fix_pdb(), axis=1)
    ligands.apply(lambda x: x.m.delete_mol(), axis=1)
    write_footer()

if __name__ == "__main__":
    main()



