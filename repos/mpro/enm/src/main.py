#!/usr/bin/env python
""" This is the master script for recreating the results.

    It imports each of the key module fucntions 
    and runs them one by one.

    Run the whole thing from the root directory 
    to replicate all the results:
    
    $ pyhton src/main.py
    or
    $ python -m src.main
"""
import get_pdb
import simulate
import clean_data
import analyze
import plot
import utilities as utils

dist_cutoff = 8.5
protein_name = "Mpro"

config = utils.read_config()
# utils.clean()

# get_pdb.main()
# simulate.main()
clean_data.main()
analyze.main(dist_cutoff)
plot.main(dist_cutoff, protein_name)

