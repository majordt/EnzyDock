#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
# To test stand alone: /home/qnt/shani/anaconda3/bin/python pdb_reader.py 4okz.pdb 4okz.pdb
# This script built mainly for TPS family, note that work fine with others ligands and cofactors. 
# 
# Initial version taken from CHARMM-GUI tutorial
# Copyright © 2022 Dan T. Major

from pymol import cmd
import sys
sys.path.append("/home/qnt/shani/anaconda3/bin/pymol")
import os
import time

class Atom:
    def __init__(self, line):
        self.lineID = line[0:6].strip()
        self.serial = int(line[6:11])
        self.name = line[11:16].strip()
        self.altLoc = line[16:17].strip()
        self.resName = line[17:21].strip()   # Include field 21 to allow 4th letter (e.g., TIP3)
        self.chainID = line[21:22]
        self.resSeq = int(line[22:26])
        self.iCode = line[26:27].strip()
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.segID = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        if self.occupancy: self.occupancy = float(self.occupancy)
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]

class Remark:
    def __init__(self, line):
        self.remark = line.strip()

class PDB:
    def __init__(self, file1):#, file2):
        self.file = file1
        self.atoms = []
        self.remarks = []
        self.parse()

    def parse(self):
        MODEL = None
        f = open(self.file, 'r')
        for line in f.readlines():
            if line.startswith('REMARK'):
                remark = Remark(line)
                self.remarks.append(remark)
            if line.startswith('MODEL'): MODEL = int(line.split()[1])
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = Atom(line)
                atom.MODEL = MODEL
                self.atoms.append(atom)
        f.close()

    def get_atoms(self, to_dict=True):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        if to_dict: return [x.__dict__ for x in self.atoms]
        else: return self.atoms

    def get_model(self, model_num, to_dict=True):
        """Return all atoms where MODEL == model_num"""
        model_atoms = [x for x in self.atoms if x.MODEL == model_num]
        if to_dict:
            return [atom.__dict__ for atom in model_atoms]
        else:
            return model_atoms

    def add_element(self):
        for atom in self.atoms:
           if atom.name in ['LI', 'BE', 'NA', 'MG', 'CA', 'SI', 'CL', 'BR', 'RB', 'CS', 'BA', 'ZN', 'CD']:
              # Special treament of carbon atoms to avoid messing up atom type in last column
              if atom.name == 'CA' or atom.name == 'CD' or atom.name == 'CL' or atom.name == 'CD':
                 if atom.resName in ['CAL', 'CAD', 'CLA', 'CAD']:
                    atom.element = atom.name[0:2]
              else:
                 atom.element = atom.name[0:2]
           else:
               atom.element = atom.name[0] # Need to define special cases for Li, Be, Na, Mg, Ca, Si, Cl, Br

    def write_pdb(self, file1):
       self.file = file1.split('.')[0] # file without .pdb 
       print(self.file)
       f_protein = self.file +'_1' +'.pdb'
       fp = open(f_protein, 'w')
       l_protein = self.file +'_ligand' + '.pdb'
       fl = open(l_protein, 'w')
       w_protein = self.file + '_wat' +'.pdb'
       fw = open(w_protein, 'w')
       line = []
       for remark in self.remarks:
           print(remark.remark, file = fp)
           print(remark.remark, file = fl)
           print(remark.remark, file = fw)
       i = 1
       j = 1
       k = 0
       resSeq = 9999
       idx = jdx = kdx = []
       for atom in self.atoms:
           if atom.lineID == 'HETATM':
               if atom.resName == 'HOH':
                   atom.resSeq = j
                   atom.resName = 'TIP3'
                   if atom.name == 'O': atom.name = 'OH2'
                   j += 1
               else:
                   if resSeq != atom.resSeq:   # New residue
                       k += 1
                       resSeq = atom.resSeq
                   atom.resSeq = k
               #k += 1
           else: # Now fix protein stuff for CHARMM
               if atom.resName == 'ILE':
                   if atom.name == 'CD1':
                       atom.name = 'CD'
               if atom.resName == 'HIS':
                   atom.resName = 'HSD'
               if atom.name == 'OXT':
                   atom.name = 'OT2'

           line = "{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"\
                  .format(atom.lineID,atom.serial,atom.name,atom.altLoc,atom.resName,atom.chainID,atom.resSeq,atom.iCode,
                          atom.x,atom.y,atom.z,atom.occupancy,atom.tempFactor,atom.element,atom.charge)

           # CHARMM style
           #line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:>2s}{:2s}"\
           #       .format(field,atom.serial,atom.name,atom.altLoc,atom.resName,atom.chainID,atom.resSeq,atom.iCode,
           #               atom.x,atom.y,atom.z,atom.occupancy,atom.tempFactor,atom.segID,atom.element,atom.charge)
           i += 1
           if atom.lineID == 'ATOM':
               print(line, file = fp)
               idx = [i, atom.resName, atom.chainID, atom.resSeq]
           elif atom.lineID == 'HETATM':
               line = line.replace('HETATM', 'ATOM  ')
               if atom.resName == 'TIP3':
                   print(line, file = fw)
                   jdx = [j - 1, atom.resName, atom.chainID, atom.resSeq]
               else:
                   #print(line)
                   print(line, file = fl)
                   kdx = [k, atom.resName, atom.chainID, atom.resSeq]
       
       field = 'TER'
       if len(idx) == 4:
          line = "{:6s}{:5d}      {:3s} {:1s}{:4d}".format(field, idx[0], idx[1], idx[2], idx[3])
          print(line, file = fp)
       if len(jdx) == 4:
          line = "{:6s}{:5d}      {:3s} {:1s}{:4d}".format(field, jdx[0], jdx[1], jdx[2], jdx[3])
          print(line, file = fw)
       if len(kdx) == 4:
          line = "{:6s}{:5d}      {:3s} {:1s}{:4d}".format(field, kdx[0], kdx[1], kdx[2], kdx[3])
          print(line, file = fl)
       line = "END"
       print(line, file = fp)
       print(line, file = fl)
       print(line, file = fw)
       fp.close()
       fl.close()
       fw.close()

    def correct_seqres(self,file1):
        self.file = file1.split('.')[0]  # file without .pdb
        print(self.file)
        f_pop = self.file + '_pop' + '.pdb'
        f_pop = open(f_pop, 'w')

        i = 1
        j = 1
        for atom in self.atoms:
            atom.resSeq = j
            line = "{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}" \
            .format(atom.lineID, atom.serial, atom.name, atom.altLoc, atom.resName, atom.chainID, atom.resSeq,
                    atom.iCode,
                    atom.x, atom.y, atom.z, atom.occupancy, atom.tempFactor, atom.element, atom.charge)
            i +=1
            if atom.lineID == 'ATOM':
                print(line, file=f_pop)
                jdx = [i, atom.resName, atom.chainID, atom.resSeq]
            field = 'TER'
        if len(jdx) == 4:
            line = "{:6s}{:5d}      {:3s} {:1s}{:4d}".format(field, jdx[0], jdx[1], jdx[2], jdx[3])
            print(line, file=f_pop)
        line = "END"
        print(line, file=f_pop)
        f_pop.close()


    def split_files(self,file1):
        protein_name = file1.split('_')[0]
        #self.file = file1.split('.')[0]
        resname_dict = {}
        for atom in self.atoms:
            #print(atom.resName)
            resname_dict[atom.chainID+atom.resName] = atom.resName
            #print(resname_dict)

        #print('key value: ', key)
        for key in resname_dict:
            fpymol = open('pymol.pml', 'w')
            print('load ' + file1, file=fpymol)

            if 'MG' not in key and 'TIP3' not in key and 'MN' not in key: #For POP- because POP is not usually in the PDB as POP
                print('select p, not name o* and not name p*', file=fpymol)
                print('remove p ', file=fpymol) #delete the C chain -if there are more need to consider here
                print('select resname ' + resname_dict[key], file=fpymol)
                print("alter (sele), resn='POP'", file=fpymol)
                print('save '+protein_name+'_'+key+'.pdb'  +','+ 'resname POP' + ' and chain '+key[0], file=fpymol)
                pdb = PDB(protein_name+'_'+key+'.pdb' )  # , sys.argv[2])
                pdb.get_atoms(to_dict=False)
                pdb.correct_seqres(protein_name+'_'+key+'.pdb')

            else: #For MG and water
                if 'MN' in key:
                    print('select resname ' + resname_dict[key], file=fpymol)
                    print("alter (sele), resn='MG'", file=fpymol)
                    print("alter (sele), elem='MG'", file=fpymol)
                    print('save '+protein_name+'_'+key+'.pdb'  +','+ 'resname MG' + ' and chain '+key[0], file=fpymol)

                else:
                    print('save '+protein_name+'_'+key+'.pdb'  +','+ 'resname ' + resname_dict[key]+ ' and chain '+key[0], file=fpymol)
            #print(key)

            fpymol.close()
            os.system('pymol -c pymol.pml')


    def split_protein(self,file1):
        protein_name = file1.split('_')[0]
        #self.file = file1.split('.')[0]
        resname_dict = {}
        for atom in self.atoms:
            #print(atom.resName)
            resname_dict[atom.chainID] = atom.chainID
            #print(resname_dict)

        #print('key value: ', key)
        for key in resname_dict:
            fpymol = open('pymol.pml', 'w')
            print('load ' + file1, file=fpymol)
            print('save '+protein_name+'_'+key+'.pdb'  +','+  ' chain '+key[0], file=fpymol)
            fpymol.close()
            os.system('pymol -c pymol.pml')

def write_header():
    print("\nPDB handle program for EnzyDock\nCopyright © 2022 Dan T. Major\n")
    print("Initial version of code taken from © CHARMM-GUI\n")

def write_footer():
    print("\nEnd PDB handle program, returning to EnzyDock main\n")

def main():
    import sys

    write_header()
    pdb = PDB(sys.argv[1])#, sys.argv[2])
    atoms = pdb.get_atoms(to_dict=False)
    pdb.add_element()
    pdb.write_pdb(sys.argv[2])

    protein = sys.argv[2].split('.')[0]

    #split ligands
    ligand_file = protein +'_ligand.pdb'
    pdb = PDB(ligand_file)  # , sys.argv[2])
    pdb.get_atoms(to_dict=False)
    pdb.split_files(ligand_file)


    #split water

    water_file = protein +'_wat.pdb'
    pdb = PDB(water_file)  # , sys.argv[2])
    pdb.get_atoms(to_dict=False)
    pdb.split_files(water_file)

    
    #split protein
    protein_file = protein +'_1.pdb'
    pdb = PDB(water_file)  # , sys.argv[2])
    pdb.get_atoms(to_dict=False)
    pdb.split_protein(protein_file)

    #for enzydock if taking just the first chain
    #dir = protein + '_final/'
    #os.system('mkdir '+ dir )
    #print(dir)
    os.system('cp '+ protein +'_ATIP3'+'.pdb' +' ' + protein +'_wat'+'.pdb' )
    os.system('cp ' + protein + '_AMG' + '.pdb' + ' ' + protein + '_mg' + '.pdb')
    os.system('cp ' + protein + '_A' + '.pdb' + ' ' + protein + '_1' + '.pdb')
    os.system('cp ' + protein + 'A*_pop' + '.pdb' + ' '+ protein + '_pop' + '.pdb') #specific here for ligand A3E9
    write_footer()

if __name__ == '__main__':
    main()
