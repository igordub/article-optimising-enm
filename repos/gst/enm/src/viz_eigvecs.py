# Plots B-factors column in PDB file as spectrum
# B-factors are replced with RMSD for each mode
#
# Usage: pymol -qc src/viz_eigvecs.py > /dev/null 

from glob import glob

cmd.delete('all')
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
    cmd.spectrum('b', 'blue_red', '%{}'.format(rmsd_name), 0.0, 1.0)

    cmd.delete(plot_name)

    cmd.group(mode_name, members=eigvecs_name)
    cmd.group(mode_name, members=rmsd_name)


cmd.set('sphere_scale', 0.6)
cmd.show_as('spheres')
cmd.show('cgo')

cmd.set_view((\
     0.781361341,   -0.511385441,    0.357710421,\
     0.429300904,    0.856461704,    0.286658257,\
    -0.452957839,   -0.070418254,    0.888746917,\
     0.000000000,    0.000000000, -160.304840088,\
    78.177932739,   63.766319275,    0.000000000,\
   126.385559082,  194.224121094,  -20.000000000 ))

cmd.set('grid_mode', 1)

cmd.bg_color('black')

cmd.save('scratch/eigvecs.pse')
