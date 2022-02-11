import pytraj as pt
from pca_mw import pca_mw

from matplotlib import pyplot as plt
import os
from os.path import join as join_paths, basename as get_basename
import numpy as np


def compute_pca(traj, enm_mask):
    """ Calculates mass-weighted PCA for trajectories.
    """
    pca_data = pca_mw(traj, mask=enm_mask,
                      fit=True, n_vecs=-1)
    _, eigdata = pca_data
    eigvals = eigdata[0]

    return eigvals


def analyze_traj(traj_path, top_path, ref_traj_path, ref_top_path, output_dir):
    """ Analyses given trajectory and compares
        it with a reference strucuture.
    """

    # Each frame is 0.01 ns
    print("Access all trajectory frames:")
    traj_all = pt.iterload(traj_path,
                           top_path,
                           frame_slice=(0, -1))
    print(traj_all)

    print("Access trajectory frames for analysis:")
    traj = pt.iterload(traj_path,
                       top_path,
                       frame_slice=(-5000, -1))
    print(traj)

    print("Reference strucutre trajectory:")
    ref_crd = pt.load(ref_traj_path, ref_top_path)
    print(ref_crd)

    backbone_residues = [2, 199, 201, 399]
    backbone_mask = "((:{}-{})|(:{}-{}))&(@CA)".format(*backbone_residues)

    enm_residues = [3, 199, 202, 398]
    enm_mask = "((:{}-{})|(:{}-{}))&(@CA)".format(*enm_residues)

    # ANALYSIS
    np_kwargs = {'fmt': ['%d','%.6e'],
                 'delimiter': ',',
                 'newline': '\n',
                 'comments':''}

    # RMSD
    rmsd_data = pt.rmsd(traj_all, mask=backbone_mask, ref=ref_crd,
        nofit=False, 
        mass=False)

    frame_nums = np.arange(traj_all.n_frames) + 1
    rmsd_data = np.vstack ((frame_nums, rmsd_data)).T

    np.savetxt(join_paths(output_dir, 'rmsd.csv'), rmsd_data,
               header="frame_number,rmsd",
               **np_kwargs)

    # RMSF
    rmsf_data = pt.rmsf(traj,
                        mask=backbone_mask)
    np.savetxt(join_paths(output_dir, 'rmsf.csv'), rmsf_data,
               header="atom_number,rmsf",
               **np_kwargs)

    # B-FACTORS
    bfactors_data = pt.bfactors(traj,
                                mask=enm_mask,
                                byres=True)

    np.savetxt(join_paths(output_dir, 'bfactors.csv'), bfactors_data,
               header="residue_number,bfactor",
               **np_kwargs)

    # EIGENVALUES
    eigvals = compute_pca(traj, enm_mask)

    modes = np.arange(len(eigvals)) + 1
    eigvals = np.vstack((modes, eigvals)).T

    np.savetxt(join_paths(output_dir, 'eigvals.csv'), eigvals,
               header="mode,eigval",
               **np_kwargs)


def main():
    """ Master script.
    """
    input_dir = 'data/05-analysis'
    output_dir = 'scratch'

    traj_path = join_paths(input_dir, 'traj.no_water.nc')
    top_path = join_paths(input_dir, 'traj.no_water.parm7')

    ref_top_path = 'data/00-structure/complex.parm7'
    ref_traj_path = 'data/00-structure/complex.ncrst'

    analyze_traj(traj_path,
                 top_path,
                 ref_traj_path,
                 ref_top_path,
                 output_dir)


if __name__ == "__main__":
    main()
