# Plots B-factors column in PDB file as spectrum
# B-factors are replced with RMSD for each mode
#
# Usage: pymol -qc src/viz_eigvecs.py > /dev/null


from pymol import cmd
from glob import glob

cmd.delete('all')
# Load modevectors module
cmd.run('src/modevectors.py')

eigvecs_filenames = glob("data/raw/project_eigvecs/mode.m????.pdb")
rmsd_filenames = glob("data/raw/project_eigvecs/rmsd.m????.pdb")

eigvecs_filenames.sort()
rmsd_filenames.sort()

if len(eigvecs_filenames) != len(rmsd_filenames):
    print("Number of RMSD files doesn't match number of eigenvector files")
    exit


for mode_idx in range(len(eigvecs_filenames)):
    eigvecs_filename = eigvecs_filenames[mode_idx]
    rmsd_filename = rmsd_filenames[mode_idx]

    mode_num = eigvecs_filename.split('/')[-1].split('.')[-2].replace('m', '')

    mode_name = 'mode_{}'.format(mode_num)
    plot_name = 'plot_{}'.format(mode_num)
    rmsd_name = 'rmsd_{}'.format(mode_num)
    eigvecs_name = 'eigvecs_{}'.format(mode_num)

    cmd.load(eigvecs_filename, plot_name)
    cmd.load(rmsd_filename, rmsd_name)

    cmd.show_as('spheres', rmsd_name)

    # Draw eigenvectors
    modevectors(rmsd_name, plot_name, outname=eigvecs_name,
                cutoff=0, cut=0, factor=1, head=0.5, head_length=2)

    # Colours residues by values in B-factor column
    cmd.spectrum('b', 'blue_white_red', '%{}'.format(rmsd_name), 0.0, 1.0)

    cmd.delete(plot_name)

    cmd.group(mode_name, members=eigvecs_name)
    cmd.group(mode_name, members=rmsd_name)


cmd.set('sphere_scale', 0.6)
cmd.show_as('spheres')
cmd.show('cgo')

origin_view = (\
     0.733638585,    0.308858454,   -0.605293334,\
    -0.673461974,    0.211624146,   -0.708282113,\
    -0.090664417,    0.927265465,    0.363261580,\
     0.000000000,    0.000000000, -189.958969116,\
   -15.172746658,  -11.892742157,   -7.302253723,\
   144.367660522,  235.550231934,   20.000000000 )

cmd.set_view(origin_view)

cmd.set('grid_mode', 1)

cmd.bg_color('black')

cmd.save('scratch/eigvecs.pse')
