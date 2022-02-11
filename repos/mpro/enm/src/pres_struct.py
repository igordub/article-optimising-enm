# PyMOL script that presents structure
#
# Usage: pymol -qc src/pres_struct.py


from pymol import cmd

cmd.delete('all')
cmd.reinitialize()

# Load complex
cmd.load('pdb/external/7bqy.pdb', 'com')

# Single letter residues
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

# Download colourblind-friendly pallete
cmd.run('https://github.com/Pymol-Scripts/Pymol-script-repo/raw/master/colorblindfriendly.py')

# Names with dots are treated special
cmd.set('group_auto_mode', 2)

# Selection
cmd.extract('prot.01', 'com and polymer and c. A')
cmd.extract('prot.02', 'com and polymer and c. B')

cmd.select('term.N.01', 'prot.01 and n. CA and i. 1')
cmd.select('term.C.01', 'prot.01 and n. CA and i. 301')
cmd.select('term.N.02', 'prot.02 and n. CA and i. 302')
cmd.select('term.C.02', 'prot.02 and n. CA and i. 602')

cmd.extract('lig.01', 'com and organic and c. A')
cmd.extract('lig.02', 'com and organic and c. B')

cmd.group('prots', 'prot.*')
cmd.group('terms', 'term.*.*')
cmd.group('ligs', 'lig.*')

cmd.delete('com prot lig')


# Representation
cmd.show_as('cartoon', 'prots')
cmd.show_as('sticks', 'ligs')
cmd.show('spheres', 'terms')


# Coloring
cmd.set('cartoon_discrete_colors', 1)

cmd.color('cb_red', 'prot.01')
cmd.color('cb_blue', 'prot.02')

cmd.color('atomic', 'ligs')
cmd.color('cb_green','ligs and elem C')


# View
origin_view = (\
    -0.088217504,    0.000000043,    0.996101260,\
    -0.000000030,    1.000000000,   -0.000000046,\
    -0.996101260,   -0.000000034,   -0.088217504,\
     0.000000000,    0.000000000, -259.065399170,\
     0.000000000,   -0.757357597,    0.000007629,\
   209.582183838,  308.548645020,  -20.000000000 )

cmd.set_view(origin_view)

cmd.center('prots')
cmd.zoom(buffer=-3,
    complete=1)

cmd.set('two_sided_lighting', 1)
cmd.set('orthoscopic', 1)
cmd.space('cmyk')

# Save session
cmd.set('ray_opaque_background', 0)

cmd.png('scratch/pres_struct.front.png', 
    width = 900,
    height = 900,
    dpi=300,
    ray=1)

cmd.turn('y',90)

cmd.png('scratch/pres_struct.side.png', 
    width = 900,
    height = 900,
    dpi=300,
    ray=1)


# Labeling
cmd.set('label_relative_mode', 0)
cmd.set('label_position', (0,0,0))

cmd.set('label_color', 'black')
cmd.set('label_outline_color', 'default')
cmd.set('label_bg_color', 'white')
cmd.set('label_bg_transparency', 0.4)

cmd.set('label_connector', 1)
cmd.set('label_connector_mode', 1)
cmd.set('label_connector_width', 4)

cmd.set('label_size', 15)


cmd.label('lig.01 and n. O1', '"N3"')
cmd.label('lig.02 and n. O1', '"N3*"')

cmd.label('term.N.01', '"N"')
cmd.label('term.C.01', '"C"')
cmd.label('term.N.02', '"N*"')
cmd.label('term.C.02', '"C*"')

# Save session
cmd.set_view(origin_view)
cmd.save('scratch/pres_struct.pse')