#!/bin/csh

if ($#argv == 0) then
   echo "Error: No arguments provided."
   echo "Usage: check_results.csh arg "
   exit
endif

# Provide paths for EnzyDock source code and target working directory
set sourcedir = #/home/qnt/majort/charmm/workspace/dock/enzydock_v1.0g6c_13092022
set workdir = #/home/qnt/majort/charmm/workspace/dock/sars-cov-2

cd $sourcedir

set sourcefiles = \
`echo crd/expl_sphere.crd \
dock/dock_ligand.inp \
dock/enzydock.inp \
dock/enzydock_main.inp \
dock/gbsw_cluster.inp \
dock/ligand_mk_conformers.inp \
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
scripts/0README \
scripts/calc_rmsd_perm.py \
scripts/cluster_ligands.py \
scripts/cluster_ligands_perm.py \
scripts/cluster_water.py \
scripts/define_system_var.csh \
scripts/files_lock.csh \
scripts/files_lock.sh \
scripts/get_qm_atom_type.py \
scripts/pdb_reader.py \
scripts/prepare.py \
scripts/python_wrapper.sh \
scripts/read_write_pdb.py \
scripts/ring_aromatic.py \
scripts/setup_enzydock_.csh \
scripts/smiles_tors.py \
scripts/softlink.py \
scripts/superpose_cluster_ligands.py \
scripts/write_qchem_inp.py \
stream/back_clear_restraints.str \
stream/back_flex.str \
stream/back_flex_ng.str \
stream/check_natoms.str \
stream/clear_all_restraints.str \
stream/clear_restraints.str \
stream/cofact_patch_offgrid.str \
stream/cofact_patch_ongrid.str \
stream/consensus_docking.str \
stream/cov_generate_offgrid.str \
stream/cov_generate_ongrid.str \
stream/cov_ligcofac.str \
stream/cov_mov_lig.str \
stream/define_flex_offgrid.str \
stream/define_flex_ongrid.str \
stream/delete_flex.str \
stream/disulf.str \
stream/expl_water.str \
stream/generate_cofactors.str \
stream/generate_grid.str \
stream/generate_protein.str \
stream/mc_tors.str \
stream/mc_trot.str \
stream/mc/mc_sidechains.str \
stream/mc/mc_sidechains_delete.str \
stream/no_samd.str \
stream/param.str \
stream/param_set_check.str \
stream/patch_offgrid.str \
stream/patch_ongrid.str \
stream/read_flex.str \
stream/read_grid.str \
stream/read_ligand.str \
stream/read_protein.str \
stream/read_water_expl.str \
stream/read_water_in.str \
stream/read_water_out.str \
stream/read_water.str \
stream/read_write_protein.str \
stream/rseed.str \
stream/samd.str \
stream/seed.str \
stream/setup_protein.str \
stream/softlink.str \
stream/top_prm_all.str \
stream/usegbsw.str \
stream/usepbeq.str \
stream/useqmmm.str \
stream/write_flex.str \
stream/zero_consensus_docking.str \
stream/consdef/consensus_def.str \
stream/consdef/userrestraints.str \
stream/gbsw/radius_gbsw.str \
stream/patching/cov_patch_g.str \
stream/patching/cov_patch_ng.str \
stream/patching/set_prot_patch.str \
stream/pbeq/radius_pbeq.str`

cd $sourcedir
echo "Current version of EnzyDock:"
pwd
ls -1 $sourcefiles
cd $workdir

foreach targetdir ( *_$argv[1] )
   cd $sourcedir
   echo "Copying all files to $workdir/$targetdir ..."
   chmod u+w -R $workdir/$targetdir
   # Use full path to override prompt
   /bin/cp --parents $sourcefiles $workdir/$targetdir
end

