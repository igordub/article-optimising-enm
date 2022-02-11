# Optimising Elastic Network Models for Protein Dynamics and Allostery: Spatial and Modal Cut-offs and Backbone Stiffness

**Authors**: [Igors Dubanevics](https://github.com/igordub) and [Tom C. B. McLeish](https://www.york.ac.uk/physics/people/mcleish/)

# Repository structure
```
├── README.md           <- The top-level README for this project.
├── ccenv.yml           <- A Conda environment file.
├── config.yaml         <- A configuration file for the scripts.
|
├── data                <- Processed data from MD and ENM simualtions.
│   ├── cap                 <- Catabolite Activator Protein data.
│   ├── gst                 <- Glutathione S-transferase data.
│   └── mpro                <- SARS-CoV-2 main protease data.
│
├── repos               <- Directorires with repository strucutres for protein MD and ENM simulations.
│   ├── cap                 <- Catabolite Activator Protein.
|   |   ├── enm                 <- ENM repository structure.
|   |   └── md                  <- MD repository structure.
│   |
│   ├── gst                 <- Glutathione S-transferase.
|   |   ├── enm
|   |   └── md
│   |
│   └── mpro                <- SARS-CoV-2 main protease.
|       ├── enm
|       └── md
|
├── fig             <- Generated graphics and figures to be used in reporting
|
├── tab             <- Processed tables to be used in reporting
|
└── src                 <- Scripts directory.
    ├── __init__.py
    ├── main.py                 <- Master script.
    ├── plot.py                 <- Vizualization script.
    └── utilities.py            <- Helpers funcitons.

```
Molecular dynamics (MD) and elastic network model (ENM) simulations are consuming and require tunning for each protein `repos/` directory contains repository structures for MD and ENM simualtions for each protein. Additionally, each of those file structures has a `README.md` file with detailed instructions on how to reproduce the simualtions.
## Re-running the example project
To re-run the plotting part of the project with conda, first install the required packages in a new conda environment:
```
conda env create -f ccenv.yml
```
You then need to activate the environment to use in your IDE. The command for this is ```conda activate comp-biophys```. Note that, to try to ensure compatibility with all operating systems, the environment file is simplified. It only specifies the version of Python, conda is trusted to resolve dependencies and pick the right versions of other required packages.

There are a couple of ways to run the python files. The first is to run them from the command line, from the root directory. 
To run all of them at once, use
``python ./main.py``
or 
``python -m src.main``
Or, from within an IDE, run the contents of src/main.py (but run it as if from the root folder).