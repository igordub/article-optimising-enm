import pytraj as pt

import os
from os.path import join as join_paths, basename as get_basename


def clean_traj(traj_paths, topology_path, output_dir):
    """ Loads production run trajectories,
        combines, strips water, autoimages,
        and saves the resulting trajectories
        and topology.
    """
    # LOAD DATA
    # Combine trajectories
    traj = pt.load(traj_paths, topology_path)
    
    print(traj)

    # CLEAN DATA
    traj = pt.strip(traj, ':WAT')
    traj.top = pt.strip(traj.top, ':WAT')

    traj = pt.autoimage(traj)
    traj = pt.superpose(traj)

    print(traj)

    # SAVE DATA
    pt.save(join_paths(output_dir, 'traj.no_water.nc'), 
            traj, overwrite=True)
    pt.save(join_paths(output_dir, 'traj.no_water.parm7'),
            traj.top, overwrite=True)

    return None

def main():
    """ Master script.
    """
    input_dir = 'data/04-production'
    traj_paths = join_paths(input_dir, 'prod.*.nc')
    topology_path = 'data/00-structure/complex.parm7'
    output_dir = 'data/05-analysis'

    clean_traj(traj_paths, topology_path, output_dir)    

if __name__ == "__main__":
    main()
    