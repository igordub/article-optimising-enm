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

cmd.set_view((\
    -0.494367570,    0.094866306,    0.864060462,\
     0.866253734,    0.136262566,    0.480662048,\
    -0.072140358,    0.986118555,   -0.149543166,\
     0.000000000,    0.000000000, -206.666610718,\
   -14.833507538,  -11.688899994,   -6.689447403,\
   162.937530518,  250.395690918,  -20.000000000 ))

cmd.set('grid_mode', 1)

cmd.bg_color('grey')

cmd.save('scratch/benms.pse')
