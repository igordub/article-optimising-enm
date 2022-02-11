import pytraj as pt
from os.path import join as join_paths, basename as get_basename

input_dir = 'data/04-production'
output_dir = 'tmp'

topology_path = 'data/00-structure/complex.parm7'
ref_crd_path = 'data/00-structure/complex.ncrst'

traj_paths = ['data/00-structure/complex.ncrst',
    'data/01-minimisation/min.01.ncrst', 'data/01-minimisation/min.02.ncrst',
    'data/02-heating/heat.nc', 'data/03-equilibration/eq.nc']

# traj_paths = join_paths('data/04-production', 'prod.0[3-9]?.nc')

traj = pt.load(traj_paths, topology_path)
print(traj)

traj = pt.strip(traj, ':WAT')
traj.topology = pt.strip(traj.topology, ':WAT')

traj = pt.autoimage(traj)
traj = pt.superpose(traj)

print(traj)
pt.save(join_paths(output_dir, 'traj.no_water.nc'), traj, overwrite=True)
pt.save(join_paths(output_dir, 'complex.no_water.parm7'),
        traj.topology, overwrite=True)
