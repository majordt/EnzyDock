#!/home/qnt/majort/anaconda3/bin/python3.6

import sys

def write_header():
    print("\nQ-Chem input writer program for EnzyDock\nCopyright © 2020 Dan T. Major\n")

def write_file1(file, qm_charge):
    file.write(
'''$comment
qchem control file
qm/mm EnzyDock simulations
$end

$rem
method                hf3c
basis                 minix
print_general_basis   true
qm_mm                 true
jobtype               sp
symmetry              off
sym_ignore            true
print_input           false
qmmm_print            true
max_scf_cycles        1000
thresh                8
scf_convergence       5
scf_guess             sad !gwh
$end

$molecule\n''')

    line = str(qm_charge) + "  1\n"
    file.write(line)
    line = "$end\n"
    file.write(line)

def write_file2(file, qm_charge):
    file.write(
'''$comment
qchem control file
EnzyDock qm/mm simulations
$end

$rem
method                hf3c
basis                 minix
print_general_basis   true
qm_mm                 true
jobtype               force
symmetry              off
sym_ignore            true
print_input           false
scf_guess             read
max_scf_cycles        1000
thresh                8
scf_convergence       5
$end

$molecule\n''')

    line = str(qm_charge) + "  1\n"
    file.write(line)
    line = "$end\n"
    file.write(line)

def main():
    write_header()
    qm_charge = sys.argv[1]
    file1 = "../stream/qmmm/qchem1.inp"
    file2 = "../stream/qmmm/qchem2.inp"
    print("Writing Q-Chem input files: \n", file1, "\n", file2, "\n")
    file = open(file1,"w")
    write_file1(file, qm_charge)
    file = open(file2,"w")
    write_file2(file, qm_charge)

if __name__ == "__main__":
    main()


