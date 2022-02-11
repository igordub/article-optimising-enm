#!/usr/bin/env python
""" This is the master script for recreating the results
    from raw data.

    It imports each of the key module fucntions 
    and runs them one by one.

    Run the whole thing from the root directory 
    to replicate all the results:
    
    $ pyhton src/05-analysis/main.py
    or
    $ python -m src.05-analysis.main
"""
import clean_traj
import analyze_traj
import plot_traj
import utilities as utils


config = utils.read_config()
# utils.clean()

# clean_traj.main()
# analyze_traj.main()
plot_traj.main()
