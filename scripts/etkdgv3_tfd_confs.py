#!/home/qnt/renanak/miniconda3/envs/my-rdkit-env/bin/python

from rdkit.Chem import rdDistGeom, AllChem, rdMolTransforms, rdMolAlign
from rdkit.Chem import Descriptors, TorsionFingerprints, rdFMCS, rdMolDescriptors
from rdkit import Chem
import rdkit
import importlib
from rdkit.Geometry import Point3D
from rdkit.Chem.Draw import rdDepictor
import numpy as np
from rdkit.ML.Cluster import Butina
import shutil
from pathlib import Path
import sys
from openbabel import pybel
import os

def find_undef_stereo_double_bonds_with_Hs(mol):
    """
    Finds double bonds in a molecule that lack defined stereochemistry (i.e., no cis/trans),
    and where this is likely due to one of the atoms in the bond being bonded only to hydrogens.
    
    For each such bond, returns:
      - The bond index
      - The index of the atom with only hydrogen neighbors
      - A list of the indices of those hydrogen atoms

    Parameters:
        mol (rdkit.Chem.rdchem.Mol): An RDKit molecule object (with explicit Hs added).
    
    Returns:
        list of tuples: Each tuple is (bond_idx, atom_idx, list_of_H_indices)
    """
    issues = []
    for bond in mol.GetBonds():
        # Only interested in double bonds with no stereochemistry
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetStereo() == Chem.BondStereo.STEREONONE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            # Check if begin atom has only hydrogens (besides its double bond partner)
            h_neighbors_begin = get_H_neighbors(begin_atom, end_atom.GetIdx())
            if len(h_neighbors_begin) == begin_atom.GetDegree() - 1:
                issues.append((bond.GetIdx(), begin_atom.GetIdx(), h_neighbors_begin))
                continue  # skip checking end_atom to avoid duplicates

            # Check if end atom has only hydrogens
            h_neighbors_end = get_H_neighbors(end_atom, begin_atom.GetIdx())
            if len(h_neighbors_end) == end_atom.GetDegree() - 1:
                issues.append((bond.GetIdx(), end_atom.GetIdx(), h_neighbors_end))

    return issues

def replace_hs_with_heavy_preserve_confs(mol, h_indices, new_atomic_num=9):
    """
    Replace given H atoms in mol with a semi-heavy atom (e.g., F) and
    preserve all conformers.

    Parameters:
        mol (rdkit.Chem.Mol): Molecule with conformers and explicit Hs.
        h_indices (list of int): Atom indices of Hs to replace.
        new_atomic_num (int): New atomic number (default: 9 for F).

    Returns:
        rdkit.Chem.Mol: New molecule with atoms replaced and conformers preserved.
    """
    rw_mol = Chem.RWMol(mol)

    # Replace H atoms
    for idx in h_indices:
        atom = rw_mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 1:
            atom.SetAtomicNum(new_atomic_num)

    new_mol = rw_mol.GetMol()

    # Copy conformers
    new_mol.RemoveAllConformers()
    for conf in mol.GetConformers():
        new_conf = Chem.Conformer(conf)
        new_mol.AddConformer(new_conf, assignId=True)

    Chem.SanitizeMol(new_mol)
    return new_mol

def get_H_neighbors(atom, exclude_idx):
    """
    Returns the indices of all hydrogen neighbors of an atom,
    excluding the atom with index `exclude_idx` (i.e., the atom
    at the other end of the double bond).
    
    Parameters:
        atom (rdkit.Chem.rdchem.Atom): The atom to check neighbors for.
        exclude_idx (int): Index of the atom to ignore (typically the bond partner).
    
    Returns:
        list of int: Indices of hydrogen neighbors.
    """
    h_indices = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() != exclude_idx and neighbor.GetAtomicNum() == 1:
            h_indices.append(neighbor.GetIdx())
    return h_indices

def cluster_conformers_by_rmsd(mol, cutoff, use_heavy_atoms=True, prealigned=False, conf_ids=None):
    """
    Clusters conformers by RMSD. Returns the grouped conformer indices
    """
    def _direct_rms_between_confs(mol, cid_i, cid_j, atom_ids=None):
        conf_i, conf_j = mol.GetConformer(cid_i), mol.GetConformer(cid_j)
        if atom_ids is None:
            atom_ids = range(mol.GetNumAtoms())
        sse = 0.0
        for a in atom_ids:
            pi, pj = conf_i.GetAtomPosition(a), conf_j.GetAtomPosition(a)
            dx, dy, dz = (pi.x - pj.x), (pi.y - pj.y), (pi.z - pj.z)
            sse += dx*dx + dy*dy + dz*dz
        return (sse / max(1, len(atom_ids))) ** 0.5

    # determine pairwise RMSD
    def _pairwise_rmsd(mol, cid_i, cid_j, atom_ids=None, prealigned=False, maxIters=50):
        if prealigned:
            return _direct_rms_between_confs(mol, cid_i, cid_j, atom_ids)
        atomMap = None if atom_ids is None else [(i, i) for i in atom_ids]
        prb = Chem.Mol(mol)  # copy so we don't touch original coords
        rms = rdMolAlign.AlignMol(prb, mol, prbCid=cid_i, refCid=cid_j, atomMap=atomMap, maxIters=maxIters)
        return float(rms)
    
    # get the lower triangular matrix of the pairwise RMSD matrix
    def _conformer_rms_utmat(mol, atom_ids=None, conf_ids=None, prealigned=False):

        if conf_ids is None:
            conf_ids = [c.GetId() for c in mol.GetConformers()]
        n = len(conf_ids)
        if n < 2:
            return [], conf_ids
        ut = []
        for i in range(n):
            for j in range(n):
                ut.append(_pairwise_rmsd(mol, conf_ids[i], conf_ids[j],
                                        atom_ids=atom_ids, prealigned=prealigned))
        N = int(np.sqrt(len(ut)))
        # Convert to condensed form by extracting the lower triangle (excluding the diagonal)
        condensed_ut = []
        for i in range(1, N):
            for j in range(i):
                condensed_ut.append(ut[i * N + j])
        assert len(condensed_ut) == n*(n-1)//2

        return condensed_ut, conf_ids
    
    if use_heavy_atoms:
        atom_ids = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    else:
        atom_ids = None
    
    ut, conf_ids = _conformer_rms_utmat(mol, atom_ids=atom_ids, conf_ids=conf_ids, prealigned=prealigned)
    n = len(conf_ids)
    if n < 2:
        return [conf_ids], ut, conf_ids
    clusters_idx = Butina.ClusterData(ut, n, cutoff, isDistData=True, reordering=False)
    clusters = [[conf_ids[k] for k in cl] for cl in clusters_idx]
    
    return clusters, ut, conf_ids

def correct_pdb(filename):
    # re-writes ligand_n.pdb with CONECT section
    for mymol in pybel.readfile("pdb", filename):
        # First backup original pdb file, then correct pdb file
        exist = True
        i = 0
        while (exist):
           bck_filename = filename + '_' + str(i)
           if (not os.path.exists(bck_filename)):
              os.system('cp -f ' + filename + ' ' + bck_filename)
              exist = False
           else:
              i += 1
        output = pybel.Outputfile("pdb", filename, overwrite=True)
        print(output.write(mymol))
        output.close()
    return bck_filename

def print_chiral_tags(mol):
    for atom in mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            print(
                f"Atom idx={atom.GetIdx()} "
                f"({atom.GetSymbol()}): {tag.name}"
            )

def calc_n_conf(mol, maxconf):
    nrot = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True)
    nconf = 3 ** nrot
    text = f"number of rotatble bonds: {nrot}"
    print(text)
    with open("enzydock.log",'a') as f:
        f.write(text + "\n")

    if maxconf > nconf:
        text = f"{maxconf} conformers requested but only {nrot} rotatble bonds.\n"
        text += f"reducing maxconf for this ligand to 3^nrot: {nconf}"
        print(text)
        with open("enzydock.log",'a') as f:
            f.write(text + "\n")
        maxconf = nconf
    return maxconf

def mk_confs(smiles, maxconf, pdb_f):
    if smiles != "":
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    else:
        # USER should make sure CONECT section is fine
        #bck_filename = correct_pdb(pdb_f)
        mol = Chem.MolFromPDBFile(
            pdb_f,
            removeHs=False,   # keep hydrogens
                # False to avoid losing info due to sanitization
        )
        Chem.AssignStereochemistryFrom3D(mol)
    print("checking maxconf")
    maxconf = calc_n_conf(mol, maxconf)
    print(Chem.MolToSmiles(mol, isomericSmiles=True))
    print("Chiral centers:", Chem.FindMolChiralCenters(mol, includeUnassigned=True))

    ps = AllChem.ETKDGv3()
    ps.numThreads = 0
    ps.trackFailures = True
    ps.useRandomCoords = True
    ps.useBasicKnowledge = True
    ps.useExpTorsionAnglePrefs = True
    #ps.randomSeed = 100
    # check all ps
    #for name in dir(ps):
    #    if name.startswith("_"):
    #        continue
    #    try:
    #        val = getattr(ps, name)
    #    except Exception:
    #        continue
    #    if callable(val):
    #        continue
    #    print(f"{name}: {val}")
    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=maxconf, params=ps)
    print(f'generated {len(cids)} conformers')

    FF='MMFF94s'  
    maxIts = 100000
    forceTol = 1e-1
    energyTol = 1e-2

    energies = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant=FF,
                                                     maxIters=maxIts)
    energies = np.array([round(e[1], 3) for e in energies])

    return energies, mol

def clust_confs(mol, energies, clusterwidth=0.2):
    tors = TorsionFingerprints.CalculateTorsionLists(mol)

    energies -= min(energies)
    ids = range(len(energies))        
    conformers = mol.GetConformers()

    # first check if there are double bonds with no stereochemistry
    # because Hs
    results = find_undef_stereo_double_bonds_with_Hs(mol)
    mod_mol = Chem.Mol(mol)
    
    # now replace those Hs with semi-heavy atom
    for bond_idx, _, h_list in results:
        mod_mol = replace_hs_with_heavy_preserve_confs(mod_mol, h_list, new_atomic_num=9)  # use F
        print(f"Modified molecule for bond {bond_idx}, Hs replaced: {h_list}")

    ## now cluster conformers to take only significant different ones
    #Cluster by cartesian rmsd if number of non-ring torsions are less than 3
    if len(tors[0]) < 3:
        print('clustering by cartesian RMSD')
        cs, _, _ = cluster_conformers_by_rmsd(mod_mol, cutoff=0.7, use_heavy_atoms=True, prealigned=False)
    else:
        dm = TorsionFingerprints.GetTFDBetweenConformers(mod_mol, ids, ids)
        
        N = int(np.sqrt(len(dm)))

        # Convert to condensed form by extracting the lower triangle (excluding the diagonal)
        condensed_dists = []
        for i in range(1, N):
            for j in range(i):
                condensed_dists.append(dm[i * N + j])
        
        cs = Butina.ClusterData(condensed_dists, len(ids), clusterwidth, 
                                isDistData=True, reordering=True)
            
    print(f'number of distinct confs: {len(cs)}')
    return cs, mol, energies

def fix_connect(to_fix_file, ref_file):
    with open(ref_file) as org_f, open(to_fix_file) as fixed_f:
        org_ls = org_f.readlines()
        fixed_ls = fixed_f.readlines()
    org_non_xyz_1 = ""
    fixed_xyz = ""
    org_non_xyz_2 = ""
    atom_tag = False
    for l in org_ls:
        if l.startswith("ATOM"):
            atom_tag = True
        elif not atom_tag:
            org_non_xyz_1 += l
        else:
            org_non_xyz_2 += l
    for l in fixed_ls:
        if l.startswith("ATOM"):
            fixed_xyz += l
    all_ls = org_non_xyz_1 + fixed_xyz + org_non_xyz_2
    with open(to_fix_file, "w") as f:
        f.write(all_ls)

#def compare_floats(a, b, precision=3):
#    return np.round(a, precision) == np.round(b, precision)
def compare_floats(a, b, precision=3):
    tol = 0.51 * 10**(-precision)
    return abs(a - b) < tol

def get_digits(xyzf, pdbf):
    xyz = {}
    pdb = {}
    crd = {}
    with open(xyzf) as f:
        xyzl = f.readlines()
    for i, line in enumerate(xyzl[2:], 1):
        atname = line[:5].strip()
        x = float(line[5:19].strip())
        y = float(line[19:33].strip())
        z = float(line[33:].strip())
        xyz[str(i)] = {"atname":atname, "x":x, "y":y, "z":z}
    with open(pdbf) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('ATOM'):
            atname = line[11:16].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            pdb[atname] = {"x":x, "y":y, "z":z}
    for pdb_atname, coors0 in pdb.items():
        for idx, coors in xyz.items():
            crd_atname = coors['atname']
            checkx = compare_floats(coors0["x"],coors["x"])
            checky = compare_floats(coors0["y"],coors["y"])
            checkz = compare_floats(coors0["z"],coors["z"])
            check_name = (pdb_atname.lower() == crd_atname.lower())
            if checkx and checky and checkz and check_name:
                crd[crd_atname] = coors
                break
    if len(crd) != len(pdb):
      print(pdbf)
    return crd

def mk_crd(crd_dict, nullcrd, out_crd):
    with open(nullcrd) as f:
        template = f.readlines()
  
    new_lines = []
    atoms = False
    for line in template:
        if atoms:
            atomno = int(line[0:10].strip())
            resno = int(line[10:20].strip())
            res = line[20:28].strip()
            type_field = line[30:38].strip()
            #x = line[40:60].strip()
            #y = line[60:80].strip()
            #z = line[80:100].strip()
            x = crd_dict[type_field]["x"]
            y = crd_dict[type_field]["y"]
            z = crd_dict[type_field]["z"]
            segid = line[102:110].strip()
            resid = line[110:118].strip()
            weighting = float(line[120:140].strip())
            formatted_line = "{:>10d}{:>10d}  {:<8s}  {:<8s}{:>20.10f}{:>20.10f}{:>20.10f}  {:<8s}  {:<8s}{:>20.10f}\n".format(
                              atomno,resno, res, type_field,x,y,z, segid, resid, weighting)
            new_lines.append(formatted_line)
        elif line.startswith("*"):
            new_lines.append(line)
        elif line[12:15] == "EXT":
            new_lines.append(line)
            atoms = True
    with open(out_crd, "w") as f:
        f.write("".join(new_lines))

def write_unique_confs(cs, mol, energies, currlig, runid, pdb_f, nullcrd):
    
    # loop over clusters (outer loop), find minimum E pose for each (inner loop)
    for i, clust in enumerate(cs):
        outname = f"../results/{currlig}_{runid}_{i+1}_0_pre"
        max_e = 1e12
        min_conf_id = clust[0]
        for conf_id in clust:
            local_e = energies[conf_id]
            if local_e < max_e:
                max_e = local_e
                min_conf_id = conf_id
        # after finding lowest e conf in cluster, write it to ED name
        # pdb
        Chem.MolToPDBFile(mol, f"{outname}.pdb", confId=min_conf_id)
        fix_connect(f"{outname}.pdb", pdb_f)
        # xyz
        conf = mol.GetConformer(min_conf_id)
        lines = []
        lines.append(str(mol.GetNumAtoms()))
        lines.append(f"conformer {i+1}")
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is not None:
                name = info.GetName().strip()
            else:
                name = atom.GetSymbol()  # fallback

            pos = conf.GetAtomPosition(atom.GetIdx())
            lines.append(f"{name:5s}{pos.x:14.8f}{pos.y:14.8f}{pos.z:14.8f}")
        
        #with open(f"{outname}.xyz", "w") as f:
        #    f.write("\n".join(lines) + "\n")
        # make sure file is fully written before continue
        tmp_xyz = f"{outname}_tmp.xyz"
        final_xyz = f"{outname}.xyz"
        
        with open(tmp_xyz, "w") as f:
            f.write("\n".join(lines) + "\n")
            f.flush()             # Python -> OS
            os.fsync(f.fileno())  # OS -> filesystem
        
        os.rename(tmp_xyz, final_xyz)  # atomic
        # crd
        crd_dict = get_digits(f"{outname}.xyz", f"{outname}.pdb")
        mk_crd(crd_dict, nullcrd, f"{outname}.crd")
        Path(f"{outname}.xyz").unlink() # delete xyz file

def fnames(currlig):
    pdb_f = f"../pdb/ligand_{currlig}.pdb"
    nullcrd = f"../crd/ligand_{currlig}.crd"
    return pdb_f, nullcrd

def write_niter(cs, currlig):
    fname = 'niters'
    charmmvar = 'maxit'
    cfilename = f"{fname}{currlig}.str"
    text = f'''* Number of clusters
*
set {charmmvar} = {len(cs)}
'''
    with open(cfilename, 'w') as f:
        f.write(text)

def handle_input():
    # Handle input arguments and return values
    # set 0 @runid
    # set 1 @currlig
    # set 2 @maxconf
    # set 3 @clusterwidth

    input_arg = sys.argv
    runid = input_arg[1].lower() # CHARMM sends in upper case parameters, change to lower
    currlig = input_arg[2].lower() # CHARMM sends in upper case parameters, change to lower
    maxconf = int(input_arg[3])
    clusterwidth = float(input_arg[4])
    #smiles = ""
    return runid, currlig, maxconf, clusterwidth #, smiles

def main():
    runid, currlig, maxconf, clusterwidth = handle_input()
    pdb_f, nullcrd = fnames(currlig)
    smiles = ""
    energies, mol = mk_confs(smiles, maxconf, pdb_f)
    cs, mol, energies = clust_confs(mol, energies, clusterwidth)
    write_unique_confs(cs, mol, energies, currlig, runid, pdb_f, nullcrd)
    write_niter(cs, currlig)

"""
def main():
    #smiles = '[CH2+]/C=C(\\C)CC/C=C(\\C)CCC=C(C)C'
    runid = "3f15"
    runid = "3f1a"
    currlig = 1
    pdb_f = f"/home/qnt/renanak/python/prepare_PL_REX/prep_ligs/pdb/010"
    pdb_f += f"/{runid.upper()}.pdb"
    crd_f = f"/home/qnt/renanak/charmm/workspace/plrex/dock/"
    crd_f += f"010_mmp12_git_aef7bde_14_12_25/{runid}_1/crd/1_{runid}_1_min_0_0.crd"
    smiles = ""
    
    #n_confs = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    n_confs = [100]#, 200, 500, 1000, 2000, 5000, 10000]
    n_clust_old = 1
    write_tag = True
    for i, n_init in enumerate(n_confs):
        cs, mol, energies = get_n_conf(smiles, n_init, pdb_f)
        n_clust = len(cs)
        if n_clust < 1.1 * n_clust_old:
            print(f"converged with {n_init} generated conformers")
            get_confs(cs, mol, energies, currlig, runid, pdb_f, crd_f)
            write_tag = False
            break
        n_clust_old = n_clust
    if write_tag:
        get_confs(cs, mol, energies, currlig, runid, pdb_f, crd_f)
"""

if __name__ == '__main__':
    main()
    #TODO: minimize in CHARMM after upload
