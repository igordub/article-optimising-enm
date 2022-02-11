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
import plot
import utilities as utils

config = utils.read_config()
# utils.clean()

plot.main()

