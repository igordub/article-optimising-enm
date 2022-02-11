""" Raw data downloader

    This script downloads a raw PDB file from PDB website.

    This script requires that packgaes within the `comp-biophys` environment 
    be installed within the Pyhton environment you are running this script in.
    See `ccenv.yml` for more deatils.

    This file can also be imported as a module and contains the following
    functions:
        * main - the main function of the script
        * run_enm - calculates an ENM
        * create_distance_matrix - calculates distance matrix for a PDB file 
        * get_cutoff_radius - finds smallest non-floppy cutoff radius for an ENM
        * brute_force_scan - performs brute-force ENM scan on PDB file

"""


import os
import utilities as utils
from urllib.request import urlretrieve

def main(output_dir):
    """ The main fucniton of this script. 

        Parameters
        ----------
        output_dir : str
            The location of the output directory

        Returns
        -------
        None

    """ 
    config = utils.read_config()
    pdb_code = config['pdb']['id']
    download_pdb(pdb_code, output_dir, biounit = True, compressed = False)

def download_pdb(pdb_code, output_dir, biounit = True, compressed = False):
    """ Downloads raw PDB files form a list of PDB IDs.
        Authored by Chris Swain (http://www.macinchem.org)
        Modified by Igors Dubanevics (https://github.com/igordub)
        Copyright CC-BY

        Parameters
        ----------
        pdb_code : str
            A four character PDB strucutre code 
        output_dir : str
            Output directory path 
        biounit : bool
            Downloads biounit. Requeired for all protomers' coordinates
            (default is True)
        compressed : bool
            Downloads compressed file (default is False) 

        Returns
        -------
        None
            The PDB file is downloaded to the ouput directory.
            
    """ 
    # Add .pdb extension and remove ':1' suffix in entities
    filename = "{:4s}.pdb".format(pdb_code[:4])
    
    # Add '1' if biounit
    if biounit:
        filename = "{}1".format(filename)
    # Add .gz extenison if compressed
    elif compressed:
        filename = "{}.gz".format(filename)
    
    url = os.path.join("https://files.rcsb.org/download/", filename.lower())
    destination_file = os.path.join(output_dir, filename)
    # Download file
    urlretrieve(url, destination_file)

    return None


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main_commandline()
