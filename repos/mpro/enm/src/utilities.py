#!/usr/bin/env python
"""
This script provides useful funcs to all other scripts
"""
import yaml
import os

def read_config():
    # Read in config file
    with open("config.yaml") as yaml_file:
        # YAML loads a list of dictionaries
        config_list = yaml.full_load(yaml_file)
        # Convert list into dict
        config_dict = {key: value for dict in config_list for key, value in dict.items()}
    return config_dict
    
def read_textfile(filepath):
    """ Read file line by line.

        Parameters
        ----------
        filepath : str
            a filepath pointing to a textfile 
            to be read

        Returns
        -------
        list
            a list containing lines form the text file
    """ 
    with open(filepath) as file: # Use file to refer to the file object
        line_list = file.read().splitlines()
    
    return line_list