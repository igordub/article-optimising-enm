#!/usr/bin/env python
"""
This script provides useful funcs to all other scripts
"""
import yaml
import glob
import os


def clean(flag="all"):
    """ Cleans project directories.
    """
    subdirs = ["raw", "interim", "processed"]
    data_filepaths = []
    pdb_filepaths = []
    for subdir in subdirs:
        data = glob.glob("data/{}/*".format(subdir))
        data_filepaths.extend(data)
        pdb = glob.glob("pdb/{}/*".format(subdir))
        pdb_filepaths.extend(pdb)

    scratch = glob.glob("scratch/*")

    if flag == "all":
        for filepath in data_filepaths:
            os.remove(filepath)
        for filepath in pdb_filepaths:
            os.remove(filepath)
        for filepath in scratch:
            os.remove(filepath)

    return None
        

def read_config():
    # Read in config file
    with open("config.yaml") as yaml_file:
        # YAML loads a list of dictionaries
        config_list = yaml.full_load(yaml_file)
        # Convert list into dict
        config_dict = {key: value for dict in config_list for key, value in dict.items()}
    return config_dict
    