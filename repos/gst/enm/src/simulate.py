""" ENM simulator
"""

import subprocess
import os
from os.path import join as join_paths, basename as get_basename, splitext
import glob
import utilities as utils


def scan_dc():
    """ Performs distsance cutoff scan. 
    """
    subprocess.run(['bash', 'src/scan_dc.sh'])

    return None


def scan_benms(dist_cutoff):
    """ Performs BENMs scan. 
    """
    subprocess.run(['bash', 'src/scan_benm.sh', '{}'.format(dist_cutoff)])

    return None


def main(dist_cutoff):
    """ Master script. 
    """
    scan_dc()
    scan_benms(dist_cutoff)

    return None


if __name__ == '__main__':
    dist_cutoff = 8.00
    main(dist_cutoff)
