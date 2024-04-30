#!/home/qnt/majort/anaconda3/envs/my-rdkit-env/bin/python3.9
# Initial version taken from CHARMM-GUI tutorial
# Copyright © 2022 Dan T. Major

class Atom:
    def __init__(self, line=None):
#        print('line=',line)
#        if line is None:
#           line = "0"*80 + "\n"   # If no file being read, initialize empty line
#           print(line)
#        print(type(line))
        self.serial = int(line[6:11])
        self.name = line[11:16].strip()
        self.altLoc = line[16:17].strip()
        self.resName = line[17:20]
        self.chainID = line[21:22]
        self.resSeq = int(line[22:26])
        self.iCode = line[26:27].strip()
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = float(line[54:60]) #.strip())
        self.tempFactor = float(line[60:66]) #.strip())
        self.segID = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        if self.occupancy: self.occupancy = float(self.occupancy)
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]

class Atom2:
    def __init__(self, serial=None, name=None, altLoc=None, resName=None, \
                       chainID=None, resSeq=None, iCode=None, \
                       x=None, y=None, z=None, \
                       occupancy=None, tempFactor=None, \
                       segID=None, element=None, charge=None):
        self.serial = -1 if serial is None else serial
        self.name = "XXXX" if name is None else name
        self.altLoc = "" if altLoc is None else altLoc
        self.resName = "XXX" if resName is None else resName
        self.chainID = "A" if chainID is None else chainID
        self.resSeq = -1 if resSeq is None else resSeq
        self.iCode = "" if iCode is None else iCode
        self.x = 0.0 if x is None else x
        self.y = 0.0 if y is None else y
        self.z = 0.0 if z is None else z
        self.occupancy = 0.0 if occupancy is None else occupancy
        self.tempFactor = 0.0 if tempFactor is None else tempFactor
        self.segID = "" if segID is None else segID
        self.element = "Z" if element is None else element
        self.charge = "" if charge is None else charge

    def __getitem__(self, key):
        return self.__dict__[key]

class Remark:
    def __init__(self, line):
        self.remark = line.strip()

class PDB:
    def __init__(self, file1=None, file2=None):
        self.atoms = []
        self.remarks = []
        if file1 is not None:
           self.file = file1
           self.parse()
        else:
           self.file = None

    def parse(self):
        MODEL = None
        f = open(self.file, 'r')
        for line in f.readlines():
            if line.startswith('REMARK'):
                remark = Remark(line)
                self.remarks.append(remark)
            if line.startswith('MODEL'): MODEL = int(line.split()[1])
            if line.startswith('ATOM'):
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

    def add_atoms(self, serial=None, name=None, altLoc=None, resName=None, \
                       chainID=None, resSeq=None, iCode=None, \
                       x=None, y=None, z=None, \
                       occupancy=None, tempFactor=None, \
                       segID=None, element=None, charge=None):
        """Add atoms (use when not reading from PDB file)"""
        atom = Atom2(serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, segID, element, charge)
        self.atoms.append(atom)

    def get_model(self, model_num, to_dict=True):
        """Return all atoms where MODEL == model_num"""
        model_atoms = [x for x in self.atoms if x.MODEL == model_num]
        if to_dict:
            return [atom.__dict__ for atom in model_atoms]
        else:
            return model_atoms

    def add_element(self):
       for atom in self.atoms:
           atom.element = atom.name[0] # Need to define special cases for Li, Be, Na, Mg, Ca, Si, Cl, Br

    def write_pdb(self, file2):
       self.file = file2
       try:
          f = open(self.file, 'w')
       except:
          print("WARNING: PDB file not found")
       line = []
       for remark in self.remarks:
           print(remark.remark, file = f)
       field = 'ATOM'
       i = 1
       for atom in self.atoms:
           if len(atom.resName) <= 3:
              line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"\
                     .format(field,atom.serial,atom.name,atom.altLoc,atom.resName,atom.chainID,atom.resSeq,atom.iCode,
                             atom.x,atom.y,atom.z,atom.occupancy,atom.tempFactor,atom.element,atom.charge)
           else:
              # Allow for TIP3
              line = "{:6s}{:5d} {:>4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"\
                     .format(field,atom.serial,atom.name,atom.altLoc,atom.resName,atom.chainID,atom.resSeq,atom.iCode,
                             atom.x,atom.y,atom.z,atom.occupancy,atom.tempFactor,atom.element,atom.charge)
           # CHARMM style
           #line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:>2s}{:2s}"\
           #       .format(field,atom.serial,atom.name,atom.altLoc,atom.resName,atom.chainID,atom.resSeq,atom.iCode,
           #               atom.x,atom.y,atom.z,atom.occupancy,atom.tempFactor,atom.segID,atom.element,atom.charge)
           print(line, file = f)
           i += 1
       field = 'TER'
       line = "{:6s}{:5d}      {:3s} {:1s}{:4d}".format(field,i,atom.resName,atom.chainID,atom.resSeq)
       print(line, file = f)
       line = "END"
       print(line, file = f)
       f.close()

#def main():
#    import sys
#    pdb = PDB(sys.argv[1])
#    atoms = pdb.get_atoms(to_dict=False)
#    pdb.add_element()
#    pdb.print_atoms()

#if __name__ == '__main__':
#    main()


