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

cmd.png('scratch/pres_eigvecs.front.png', \
    width = 900,\
    height = 2700,\
    dpi=300,\
    ray=1)

cmd.turn('y',90)

cmd.png('scratch/pres_eigvecs.side.png', \
    width = 900,\
    height = 2700,\
    dpi=300,\
    ray=1)


# Show a problematic region
problem_view = (\
    -0.996026933,    0.000000068,   -0.089051567,\
     0.000000070,    0.999999940,   -0.000000016,\
     0.089051567,   -0.000000023,   -0.996026933,\
     0.000000000,    0.000000000,  -44.626552582,\
    -3.722499847,   -9.086000443,  -28.888999939,\
    18.057998657,   71.195106506,  -20.000000000 )

cmd.set_view(problem_view)

# Selection
cmd.select('prob_nodes', 'resi 492 + resi 493')

# Label
cmd.set('label_size', 50)
cmd.set('label_position', (0, 0, 2.5))
cmd.set('label_color', 'white')
cmd.set('label_shadow_mode', 3)
cmd.set('label_font_id', 10)

cmd.label('prob_nodes', '"*"')

cmd.png('scratch/pres_eigvecs.region.png', \
    width = 900,\
    height = 2700,\
    dpi=300,\
    ray=1)

# Save session
cmd.set_view(origin_view)
cmd.save('scratch/pres_eigvecs.pse')