# PyMOL script that presents eigenvectors
# problematic regions for ENM physicality
#
# Usage: pymol -qc src/pres_benms.py


from pymol import cmd


# View
cmd.set('grid_mode', 1)

origin_view = (\
     0.733638585,    0.308858454,   -0.605293334,\
    -0.673461974,    0.211624146,   -0.708282113,\
    -0.090664417,    0.927265465,    0.363261580,\
     0.000000000,    0.000000000, -189.958969116,\
   -15.172746658,  -11.892742157,   -7.302253723,\
   144.367660522,  235.550231934,   20.000000000 )

cmd.set_view(origin_view)

cmd.set('two_sided_lighting', 1)
cmd.set('orthoscopic', 1)
cmd.space('cmyk')

# Save images
cmd.set('ray_opaque_background', 0)

cmd.png('scratch/pres_benms.png', \
    width = 2700,\
    height = 900,\
    dpi=300,\
    ray=1)

# Save session
cmd.set_view(origin_view)
cmd.save('scratch/pres_benms.pse')