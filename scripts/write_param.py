import os

# extracts each word from current path
split_dir = os.getcwd().split('/')

# inserts final dir into proj and one level up path into DIR0
paramkeys = {
    'proj': split_dir[-2],
    'DIR0': '/'.join(split_dir[:-2]),
    'user': split_dir[3]
}

# prints proj and DIR values
print(f"\nThese values are being written to param.str:\nset proj {paramkeys['proj']}\nset DIR0 {paramkeys['DIR0']}\n")

# param.str content
paramstr = ['*********************************************************\n',
    '*                          param.str\n',
    '*********************************************************\n',
    '* This script streams in all important global variables\n',
    '* needed for EnzyDock.\n',
    '* For variables often changed by user, see userparam.str \n',
    "* Project name 'proj' must be set and should equal folder name.\n",
    '*\n',
    '\n',
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n',
    '! General parameters\n',
    '\n',
    f'set proj          {paramkeys["proj"]}                                ! Project name\n',
    '\n',
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n',
    '! Directories\n',
    f'set DIR0          "{paramkeys["DIR0"]}"      ! If mixed case path name, must be in quotation marks\n',
    'set softlink true\n',
    'stream ../stream/softlink.str 1\n',
    'set topDIR       @DIR/local_top                      ! Local parameter directory\n',
    'set potDIR       @DIR/local_pot                      ! Local potentials directory\n',
    'set strDIR       @DIR/stream                         ! Local stream directory\n',
    'set stpDIR       @strDIR/patching                    ! Local directory containing all patching calls\n',
    'set stdDIR       @strDIR/consdef                     ! Local directory containing specific constraints definitions\n',
    'set stmDIR       @strDIR/mc                          ! Local directory containing MC moves to add/delete\n',
    'set PSFDIR       @DIR/psf                            ! Local psf directory\n',
    'set pdbDIR       @DIR/pdb                            ! Local pdb directory\n',
    'set crdDIR       @DIR/crd                            ! Local crd directory\n',
    'set resDIR       @DIR/results                        ! Results directory\n',
    'set scrDIR       @DIR/scr                            ! Scratch disk for current project\n',
    f'set scr1DIR      /var/spool/scratch/{paramkeys["user"]}/enzydock/scr/@proj      ! Multiple project scratch\n',
    f'set scr2DIR      /var/spool/scratch/{paramkeys["user"]}/enzydock/scr/@proj      ! Multiple project scratch\n',
    'set docDIR       @DIR/dock                           ! The dock directory that holds all files\n',
    'set cgenff       /private/chem/silcsbio.2023.1/cgenff/cgenff      \n',
    '\n',
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n',
    '! RTF and PRM files\n',
    '!set ProtTopFile36   @topDIR/top_all36_prot.rtf                ! Protein topology, C36 FF (August 2012)\n',
    '!set ProtParFile36   @topDIR/par_all36m_prot.prm               ! Protein parameters, C36 FF (August 2012)\n',
    '!set WatIonTopFile   @topDIR/top_all36_water_ions_dd.rtf       ! Water and ion topology\n',
    '!set WatIonParFile   @topDIR/par_all36_water_ions_dd.prm       ! Water and ion parameters\n',
    '!set Tip3File        @crdDIR/tip216.crd                        ! Equilibrated TIP3 file\n',
    'set gridrtf         @topDIR/grid_probes.rtf                   ! rtf for grid generation on receptor\n',
    'set gridprm         @topDIR/grid_probes.prm                   ! prm for grid generation on receptor\n',
    '\n',
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n',
    '!  Initialization\n',
    '\n',
    '! Initialize random number engine\n',
    'stream @strDIR/rseed.str\n',
    '\n',
    '! Initialize all run parameters\n',
    'stream @strDIR/userparam.str\n',
    '\n',
    '! Check and initialize undefined parameters and check consistencies\n',
    'stream @strDIR/param_set_check.str\n',
    '\n',
    'return\n',
    '\n'
]

# write all lines into param.str file
with open('../stream/param.str', 'w') as f:
    f.writelines(paramstr)
