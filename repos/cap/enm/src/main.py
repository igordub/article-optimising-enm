#!/usr/bin/env python
""" This is the master script for recreating the results.

    It imports each of the key module fucntions 
    and runs them one by one.

    Run the whole thing from the root directory 
    to replicate all the results:
    
    $ python src/main.py
    or
    $ python -m src.main
"""
import get_pdb
import simulate
import clean_data
import analyze
import plot
import utilities as utils

dist_cutoff = 8.0 
protein_name = "CAP"

config = utils.read_config()
# utils.clean()

# get_pdb.main()
# simulate.main(dist_cutoff)
# clean_data.main()

# clean_data.clean_dc()
# analyze.analyze_dc(dist_cutoff)
# plot.plot_dc(dist_cutoff, protein_name)

clean_data.clean_benms()
analyze.analyze_benms(dist_cutoff)
plot.plot_benms(dist_cutoff, protein_name)