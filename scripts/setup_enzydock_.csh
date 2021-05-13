#! /bin/tcsh
# Setup enzydock run
# Should be run from project directory
# Usage: ./enzydock_setup.tcsh
# Written by DTM 12/09/2019
# Updated by DTM 26/04/2020

if ($#argv > 0) then
   echo "Error: No arguments needed."
   echo "Usage: ./enzydock_setup.tcsh"
endif

# DIR points to clean, system un-specific code
set DIR = /home/qnt/majort/charmm/workspace/dock/enzydock_v1.0e
# SCRDIR points to location where EnzyDock will write grid files
set SCRDIR = /home/qnt/majort/charmm/workspace/dock/scr

if (! -e $DIR) then
   echo "No EnzyDock directory found:"
   echo $DIR
   echo "Check DIR definition in this script and location of EnzyDock code."
   exit 1
endif

if (! -e $SCRDIR) then
   echo "No scratch directory found, creating scratch:"
   echo $SCRDIR
   mkdir $SCRDIR
endif

# Cleaning up previous EnzyDock run
if ( -e crd || -e dock || -e docs || -e local_pot || -e local_top || -e pdb || -e psf || -e results || -e scr || -e stream || -e scripts ) then
   echo "Existing local EnzyDock directory found in current directory."
   # Move system specific files out before removing
#   /bin/cp -f stream/*param.str . > & /dev/null
#   /bin/cp -f pdb/*_[1-9]*.pdb . > & /dev/null
#   /bin/cp -f local_top/*_[1-9]*.str . > & /dev/null
#   /bin/cp -f local_top/*_[1-9]*.rtf . > & /dev/null
#   /bin/cp -f local_top/*_[1-9]*.prm . > & /dev/null
#   /bin/cp -f stream/*/* . > & /dev/null
   set bool = 1
   while ($bool)
      echo "Do you want to remove old EnzyDock directory? [y/n]: "
      set yn = $<
      switch ($yn)
         case [Yy]* :
            rm -rf crd dock docs local_pot local_top pdb psf results scr stream scripts
            set bool = 0
            breaksw
         case [Nn]* :
            set bool = 0
            breaksw
         case * :
            echo "Please answer yes or no.";;
            breaksw
      endsw
   end
endif

# Copy EnzyDock version to current working directory and move user created files into correct folders
cp -r $DIR/* .
chmod u+rwx -R *
/bin/mv -f consensus_def.str         stream/consdef/  > & /dev/null
/bin/mv -f userrestraints.str        stream/consdef/  > & /dev/null
/bin/mv -f mc_add.str                stream/mc/       > & /dev/null
/bin/mv -f mc_delete.str             stream/mc/       > & /dev/null
/bin/mv -f patch_g_*.str             stream/patching/ > & /dev/null
/bin/mv -f patch_ng_*.str            stream/patching/ > & /dev/null
/bin/mv -f cov_patch_g.str           stream/patching/ > & /dev/null
/bin/mv -f cov_patch_ng.str          stream/patching/ > & /dev/null
/bin/mv -f cofact*_patch_g_*.str     stream/patching/ > & /dev/null
/bin/mv -f cofact*_patch_ng_*.str    stream/patching/ > & /dev/null
/bin/mv -f *param.str                stream/          > & /dev/null
/bin/mv -f *.pdb                     pdb/             > & /dev/null
/bin/mv -f *.str *.rtf *.prm         local_top/       > & /dev/null


