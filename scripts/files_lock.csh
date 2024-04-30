#!/bin/csh

# Write-protect all main EnzyDock files, so that they can't be changed during run
chmod -w ../crd/expl_sphere.crd \
dock_ligand.inp \
enzydock.inp \
enzydock_main.inp \
gbsw_cluster.inp \
ligand_mk_conformers.inp \
ligand_simulate.inp \
mm_cluster.inp \
mm.inp \
pbeq_cluster.inp \
qmmm_cluster.inp \
qmmm.inp \
rmsd_dock.inp \
../local_top/covalent.str \
../local_top/grid_probes.prm \
../local_top/grid_probes.rtf \
../local_top/par_all36_cgenff.prm \
../local_top/par_all36_prot.prm \
../local_top/top_all36_cgenff.rtf \
../local_top/top_all36_prot.rtf \
../local_top/toppar_water_ions.str \
../scripts/calc_rmsd_perm.py \
../scripts/cluster_ligands.py \
../scripts/cluster_ligands_perm.py \
../scripts/cluster_water.py \
../scripts/define_system_var.csh \
../scripts/files_lock.csh \
../scripts/files_lock.sh \
../scripts/get_qm_atom_type.py \
../scripts/pdb_reader.py \
../scripts/prepare.py \
../scripts/python_wrapper.sh \
../scripts/read_write_pdb.py \
../scripts/ring_aromatic.py \
../scripts/setup_enzydock_.csh \
../scripts/smiles_tors.py \
../scripts/softlink.py \
../scripts/superpose_cluster_ligands.py \
../scripts/write_qchem_inp.py \
../stream/back_clear_restraints.str \
../stream/back_flex.str \
../stream/back_flex_ng.str \
../stream/check_natoms.str \
../stream/clear_all_restraints.str \
../stream/clear_restraints.str \
../stream/cofact_patch_offgrid.str \
../stream/cofact_patch_ongrid.str \
../stream/consensus_docking.str \
../stream/cov_generate_offgrid.str \
../stream/cov_generate_ongrid.str \
../stream/cov_ligcofac.str \
../stream/cov_mov_lig.str \
../stream/define_flex_offgrid.str \
../stream/define_flex_ongrid.str \
../stream/delete_flex.str \
../stream/disulf.str \
../stream/expl_water.str \
../stream/generate_cofactors.str \
../stream/generate_grid.str \
../stream/generate_protein.str \
../stream/mc_tors.str \
../stream/mc_trot.str \
../stream/mc/mc_sidechains.str \
../stream/mc/mc_sidechains_delete.str \
../stream/no_samd.str \
../stream/param.str \
../stream/param_set_check.str \
../stream/patch_offgrid.str \
../stream/patch_ongrid.str \
../stream/read_flex.str \
../stream/read_grid.str \
../stream/read_ligand.str \
../stream/read_protein.str \
../stream/read_water_expl.str \
../stream/read_water_in.str \
../stream/read_water_out.str \
../stream/read_water.str \
../stream/read_write_protein.str \
../stream/rseed.str \
../stream/samd.str \
../stream/seed.str \
../stream/setup_protein.str \
../stream/softlink.str \
../stream/top_prm_all.str \
../stream/usegbsw.str \
../stream/usepbeq.str \
../stream/useqmmm.str \
../stream/write_flex.str \
../stream/zero_consensus_docking.str \
../stream/consdef/consensus_def.str \
../stream/consdef/userrestraints.str \
../stream/gbsw/radius_gbsw.str \
../stream/patching/cov_patch_g.str \
../stream/patching/cov_patch_ng.str \
../stream/patching/set_prot_patch.str \
../stream/pbeq/radius_pbeq.str

