# PyMOL script that presents eigenvectors
# problematic regions for ENM physicality
#
# Usage: pymol -qc src/pres_benms.py


from pymol import cmd


# View
cmd.set('grid_mode', 1)

origin_view = (
    0.872755289,    0.487586439,   -0.023609711,
    0.486824632,   -0.865779936,    0.115871824,
    0.036053561,   -0.112623222,   -0.992978692,
    0.000000000,    0.000000000, -196.234619141,
    78.121902466,   63.862365723,    0.000000000,
    157.996505737,  234.472732544,  -20.000000000)

cmd.set_view(origin_view)

cmd.turn('x',-90)

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