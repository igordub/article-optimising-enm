# PyMOL script that presents eigenvectors
# problematic regions for ENM physicality
#
# Usage: pymol -qc src/pres_benms.py


from pymol import cmd


# View
cmd.set('grid_mode', 1)

origin_view = (\
     0.089051612,    0.000000068,   -0.996026933,\
     0.000000016,    1.000000000,    0.000000070,\
     0.996026933,   -0.000000023,    0.089051612,\
     0.000000000,    0.000000000, -264.891601562,\
     0.000000000,   -0.639472961,    0.000003815,\
   219.206375122,  310.576721191,  -20.000000000 )

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