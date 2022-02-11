# Draws EN for different distance cutoffs
#
# Usage: pymol -qc src/viz_dc_enms.py > /dev/null

from glob import glob

cmd.delete('all')

caonly_path = "data/raw/scan-dc/c08.00/0/CAonly.pdb"
filenames = glob("data/raw/scan-dc/c??.??/0/draw_enm.pml")
filenames.sort()

for filename in filenames[6:14:2]:
    print(filename)
    cutoff = filename.split('/')[-3]

    cmd.load(caonly_path, 'CAonly')
    cmd.run(filename)

    cmd.set_name('CAonly', cutoff)

cmd.set_view((\
     0.733638585,    0.308858454,   -0.605293334,\
    -0.673461974,    0.211624146,   -0.708282113,\
    -0.090664417,    0.927265465,    0.363261580,\
     0.000000000,    0.000000000, -189.958969116,\
   -15.172746658,  -11.892742157,   -7.302253723,\
   144.367660522,  235.550231934,   20.000000000 ))

cmd.set('grid_mode', 1)

cmd.bg_color('grey')

cmd.save('scratch/dc_enms.pse')
