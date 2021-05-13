#!/bin/csh

set sourcedir = /home/qnt/majort/charmm/workspace/dock/enzydock_v1.0g6
set workdir = /home/qnt/majort/charmm/workspace/dock/test

cd $sourcedir

set sourcefiles = \
`echo dock/dock_ligand.inp \
dock/enzydock.inp \
dock/enzydock_main.inp \
dock/gbsw_cluster.inp \
dock/ligand_simulate.inp \
dock/mm_cluster.inp \
dock/mm.inp \
dock/pbeq_cluster.inp \
dock/qmmm_cluster.inp \
dock/qmmm.inp \
dock/rmsd_dock.inp \
local_top/covalent.str \
local_top/grid_probes.prm \
local_top/grid_probes.rtf \
local_top/par_all36_cgenff.prm \
local_top/par_all36_prot.prm \
local_top/top_all36_cgenff.rtf \
local_top/top_all36_prot.rtf \
local_top/toppar_water_ions.str \
scripts/files_lock.csh \
scripts/get_qm_atom_type.py \
scripts/calc_rmsd_perm.py \
scripts/cluster_ligands.py \
scripts/cluster_ligands_perm.py \
scripts/read_write_pdb.py \
scripts/setup_enzydock_.csh \
scripts/smiles_tors.py \
scripts/softlink.py \
scripts/write_qchem_inp.py \
stream/back_clear_restraints.str \
stream/back_flex.str \
stream/clear_restraints.str \
stream/cofact_patch_offgrid.str \
stream/cofact_patch_ongrid.str \
stream/consensus_docking.str \
stream/cov_generate_offgrid.str \
stream/cov_generate_ongrid.str \
stream/cov_mov_lig.str \
stream/define_flex_offgrid.str \
stream/define_flex_ongrid.str \
stream/delete_flex.str \
stream/disulf.str \
stream/generate_cofactors.str \
stream/generate_grid.str \
stream/generate_protein.str \
stream/mc_tors.str \
stream/mc_trot.str \
stream/mc/mc_sidechains.str \
stream/mc/mc_sidechains_delete.str \
stream/setup_protein.str \
stream/param_set_check.str \
stream/patch_offgrid.str \
stream/patch_ongrid.str \
stream/read_flex.str \
stream/read_grid.str \
stream/read_ligand.str \
stream/read_protein.str \
stream/read_water.str \
stream/read_write_protein.str \
stream/rseed.str \
stream/samd.str \
stream/seed.str \
stream/softlink.str \
stream/top_prm_all.str \
stream/usegbsw.str \
stream/usepbeq.str \
stream/useqmmm.str \
stream/write_flex.str \
stream/consdef/consensus_def.str \
stream/consdef/userrestraints.str \
stream/gbsw/radius_gbsw.str \
stream/patching/cov_patch_g.str \
stream/patching/cov_patch_ng.str \
stream/pbeq/radius_pbeq.str \`

cd $sourcedir
echo "Current version of EnzyDock:"
pwd
ls -1 $sourcefiles
cd $workdir

cd $sourcedir
echo "Copying all files to $workdir ..."
chmod u+w -R $workdir
# Use full path to override prompt
/bin/cp --parents $sourcefiles $workdir

