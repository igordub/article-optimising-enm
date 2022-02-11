# PyMOL script that presents eigenvectors
# problematic regions for ENM physicality
#
# Usage: pymol -qc src/pres_benms.py


from pymol import cmd

cmd.delete('all')
cmd.reinitialize()

# Load complex
cmd.load('scratch/benms.pse')

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
cmd.disable('not b0001 + b0010 + b0100')

# Representation
cmd.hide('spheres')
cmd.show('spheres', 'b0001')

# Coloring
# Color into tab20 second color form matplotlib
# colormaps which corresponds to 8.0 A distance cutoff
cmd.color('0xaec7e8', 'b0001 and elem C') 
