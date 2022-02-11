from pytraj import matrix
from pytraj.all_actions import mean_structure, projection, get_reference


def pca_mw(traj,
           mask,
           n_vecs=2,
           fit=True,
           ref=None,
           ref_mask=None,
           dtype='ndarray',
           top=None):
    ref_mask_ = ref_mask if ref_mask is not None else mask
    atom_indices = traj.top.select(mask)
    mass = traj.top.mass[atom_indices]

    if fit:
        if ref is None:
            traj.superpose(ref=0, mask=ref_mask_)
            avg = mean_structure(traj)
            traj.superpose(ref=avg, mask=ref_mask_)
            n_refs = 2
        else:
            ref = get_reference(traj, ref)
            traj.superpose(ref=ref, mask=ref_mask_)
            n_refs = 1

    avg2 = mean_structure(traj, mask=mask)
    mat = matrix.mwcovar(traj, mask)
    if n_vecs < 0:
        # Exclude trivial modes
        n_vecs = mat.shape[0] - 6
    else:
        n_vecs = n_vecs

    eigenvectors, eigenvalues = matrix.diagonalize(
        mat, n_vecs=n_vecs, dtype='tuple',
        scalar_type='mwcovar',
        mass=mass)
    projection_data = projection(
        traj,
        mask=mask,
        average_coords=avg2.xyz,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        scalar_type='mwcovar',
        dtype=dtype)

    # release added transformed commands for TrajectoryIterator
    if fit and hasattr(traj, '_transform_commands'):
        for _ in range(n_refs):
            traj._transform_commands.pop()
        if traj._transform_commands:
            traj._reset_transformation()
        else:
            traj._remove_transformations()
    return projection_data, (eigenvalues, eigenvectors)
