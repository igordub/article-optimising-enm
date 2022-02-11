# PyMOL script that presents eigenvectors
# problematic regions for ENM physicality
#
# Currently, grid mode is not working in batch mode
# i.e. in CLI. Use workaround, by running the script
# wiht GUI and then form PyMOL CLI save the session.
#
# Usage:    $ pymol src/pres_eigvecs.py
#           > run src/save.pres_eigvecs.py


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

cmd.set('two_sided_lighting', 1)
cmd.set('orthoscopic', 1)
cmd.space('cmyk')

# Save images
cmd.set('ray_opaque_background', 0)

cmd.png('scratch/pres_eigvecs.front.png',
        width=900,
        height=2700,
        dpi=300,
        ray=1)

cmd.turn('x', -90)

cmd.png('scratch/pres_eigvecs.bottom.png',
        width=900,
        height=2700,
        dpi=300,
        ray=1)


# Show a problematic region
problem_view = (\
     0.872755229,   -0.023609731,   -0.487586439,\
     0.486824602,    0.115871862,    0.865779936,\
     0.036053557,   -0.992978692,    0.112623267,\
     0.000014203,    0.000107004,  -42.171779633,\
    52.411293030,   62.070232391,   -5.873627663,\
     6.529368877,   77.810188293,  -20.000000000 )

cmd.set_view(problem_view)

# Selection
cmd.select('prob_nodes', 'resi 27 + resi 28')

# Label
cmd.set('label_size', 50)
cmd.set('label_position', (0, 0, 2.5))
cmd.set('label_color', 'white')
cmd.set('label_shadow_mode', 3)
cmd.set('label_font_id', 10)

cmd.label('prob_nodes', '"*"')

cmd.png('scratch/pres_eigvecs.region.png',
        width=900,
        height=2700,
        dpi=300,
        ray=1)

# Save session
cmd.set_view(origin_view)
cmd.save('scratch/pres_eigvecs.pse')
