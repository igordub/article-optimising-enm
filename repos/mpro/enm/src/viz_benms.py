# Draws EN for different backbone stiffening
# coefficients
#
# Usage: pymol -qc src/viz_benms.py > /dev/null

from glob import glob

cmd.delete('all')

caonly_path = "data/raw/scan-benm/b0001/0/CAonly.pdb"
filenames = glob("data/raw/scan-benm/b????/0/draw_enm.pml")
filenames.sort()

for filename in filenames:
    print(filename)
    back_coef = filename.split('/')[-3]

    cmd.load(caonly_path, 'CAonly')
    cmd.run(filename)

    cmd.set_name('CAonly', back_coef)

origin_view = (\
     0.089051612,    0.000000068,   -0.996026933,\
     0.000000016,    1.000000000,    0.000000070,\
     0.996026933,   -0.000000023,    0.089051612,\
     0.000000000,    0.000000000, -264.891601562,\
     0.000000000,   -0.639472961,    0.000003815,\
   219.206375122,  310.576721191,  -20.000000000 )

cmd.set_view(origin_view)

cmd.set('grid_mode', 1)

cmd.bg_color('grey')

cmd.save('scratch/benms.pse')
