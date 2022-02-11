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
     0.733638585,    0.308858454,   -0.605293334,\
    -0.673461974,    0.211624146,   -0.708282113,\
    -0.090664417,    0.927265465,    0.363261580,\
     0.000000000,    0.000000000, -202.689926147,\
   -15.172746658,  -11.892742157,   -7.302253723,\
   157.098617554,  248.281188965,   20.000000000 )

cmd.set_view(origin_view)

cmd.set('two_sided_lighting', 1)
cmd.set('orthoscopic', 1)
cmd.space('cmyk')

# Save images
cmd.set('ray_opaque_background', 0)

cmd.hide('labels')

cmd.png('scratch/pres_eigvecs.front.png', \
    width = 900,\
    height = 2700,\
    dpi=300,\
    ray=1)

cmd.turn('y',-45)

cmd.png('scratch/pres_eigvecs.side.png', \
    width = 900,\
    height = 2700,\
    dpi=300,\
    ray=1)


# Show a problematic region
problem_view = (\
     0.946767807,    0.308858454,    0.090753794,\
     0.024621546,    0.211624146,   -0.977040589,\
    -0.320974141,    0.927265465,    0.192755312,\
     0.000000000,    0.000000000,  -46.418090820,\
   -11.044500351,  -31.895999908,   12.901750565,\
    25.045850754,   67.790298462,   20.000000000 )

cmd.set_view(problem_view)

# Selection
cmd.select('prob_nodes', 'i. 396-397')

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