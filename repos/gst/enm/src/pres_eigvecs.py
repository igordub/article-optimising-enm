# PyMOL script that presents eigenvectors
# problematic regions for ENM physicality
#
# Currently, grid mode is not working in batch mode
# i.e. in CLI. Use workaround, by running the script
# wiht GUI and then form PyMOL CLI save the session.
#  
# Usage:
#           $ pymol src/pres_eigvecs.py
#      PyMOL> run src/save.pres_eigvecs.py


from pymol import cmd

cmd.delete('all')
cmd.reinitialize()

# Load complex
cmd.load('scratch/eigvecs.pse')

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
pres_modes = [9,10,11]
for mode_num in range(50):
    mode_num += 1
    if not mode_num in pres_modes:
        cmd.delete('*_{:04d}'.format(mode_num))

