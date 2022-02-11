#! /bin/bash

# Extracts slice of pfENM bonds from PyMOL script

CONDITION="\$2 > 0.125"

awk -F ", " '($1 == "cmd.set_bond(\"stick_radius\"") && ((0.04 < $2) && ($2 < 0.6)) {print f} {f=$0}' draw_enm.pml > bond.pml
awk -F ", " '($1 == "cmd.set_bond(\"stick_radius\"") && ((0.04 < $2) && ($2 < 0.6))' draw_enm.pml > set_bond.pml

grep -v 'cmd.set_bond' draw_enm.pml | grep -v 'cmd.bond' > draw_enm.slice.pml
# echo 'cmd.set_bond("stick_radius", 0.0, "all", "all")' >> draw_enm.slice.pml

cat bond.pml >> draw_enm.slice.pml
cat set_bond.pml >> draw_enm.slice.pml
