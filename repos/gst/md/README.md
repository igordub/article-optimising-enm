# All-atom Molecular Dynamics Simulation of the Glutathione S-transferase

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
    <li><a href="#references">References</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->

## About The Project

We perform all-atom MD simualtion of the Glutathione S-transferase (GST) to inform construction of the Elastic Network Model (ENM) for the same protein.

PDB ID: [1M9A](https://www.rcsb.org/structure/1M9A)

### Built With

- [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science)
- [AMEBR](https://ambermd.org)
- [AmberTools](https://ambermd.org/AmberTools.php)
- [PyTraj](https://amber-md.github.io/pytraj/latest/index.html)

<!-- GETTING STARTED -->

## Getting Started

Simply, clone the repo to your machine.
Clone the repo:

```bash
git clone https://github.com/igordub/md-gst.git
```

### Prerequisites

The scripts were tested on [the Viking](https://www.york.ac.uk/it-services/services/viking-computing-cluster/) (research computing clsuter at the University of York). Any other machine that uses [Slurm Workload Manager](https://slurm.schedmd.com/) will be suitable for the job.

#### Complex Prepartaion

AmberTools was used to prepare the strucure for simualtion.

#### Simulation

Strucutre preparation, minimisation, heating and initial equilibration can be done on a desktop using `pmemd` program form `AMBER` in parallel with MPI.
However, further simualtions require GPU computation using `pmemd.cuda` which is only available with the official license, either industrial or academic.
Ask your local IT support team if not sure.

### Installation

If you have Conda installed on your machine simply create 'amber' Conda environment from `ccenv.yml` file.

```bash
conda env create -f ccenv.yml
```

or intall quietly in the background

```bash
conda env create -f ccenv.yml > /dev/null &
```

<!-- USAGE EXAMPLES -->

## Usage

Give execution permission to all shell scritps:

```bash
chmod u+x src/**/*.sh
```

### Complex Prepartaion

1. Download the strucutre:

   ```bash
   ./src/00-structure/01.download_pdb.sh 1M9A
   ```

2. Prepare protein and ligand strucutres:

   ```bash
   ./src/00-structure/02.prepare_pdb.sh
   ```

3. Create force fileds for ligands:
   ```bash
   ./src/00-structure/02.configure_ligand.sh
   ```
4. Create coordinates and topology for the complex:
   ```bash
   tleap -sf src/00-structure/tleap.04.leapin
   ```
   If you want to check how the coordinates and topology are built step-by-step, run `src/00-structure/tleap.0[1-4].leapin` scripts consequently.
   **Note** Outputs form the first three scripts are saved in `tmp/` directory.

Now, we are ready to minimize our complex solvated in water.

### Minimization, Heating, Equilibration and Production
Each simulation stage can be performed on CPU and GPU with `job.cpu.sh` and `job.gpu.sh` `SLURM` job scripts, respectively.  
Use `sbatch` command to sumbit job scripts to `SLRUM`.

```bash
sbatch src/01-minimisation/job.gpu.sh   # Run time: 00:05:00
sbatch src/02-heating/job.gpu.sh        # Run time: 00:05:00
sbatch src/03-equilibration/job.gpu.sh  # Run time: 00:05:00
```
The equlibration script prepares system for the actual equilbiration which is done by the prodcution run scripts.
Except at least 100 ns 
Equilibraton/produciton will take a lot of computational resource. We will are aiming for 400 ns production run.
`.mdin` is set for a 5 ns run which completes in ~45 minutes on NVIDIA Tesla V100 32GB SXM2.

```bash
sbatch src/04-production/job.gpu.sh   # Run time: 00:45:00
```

You can change number of 5 ns runs in `job.{cpu,gpu}.sh` in `data/04-production` by changing commadline arguments fed into `run_production.sh`.
```bash
...
./src/04-production/run_production.sh "1" "10"
...
```
The commnad above will run will run production for 50 ns (10 runs for 5 ns each) staring with trajectory right after initial equilibration.
Change requested time in the job script accordingly to number of production runs.

## Miscellaneous

- AMBER input files **must** end with an empty line. Otherwise, Fortran executable produce the following segmentation error:

  ```
  At line 639 of file mdin_ctrl_dat.F90 (unit = 5, file = 'src/04-production/prod.mdin')
  Fortran runtime error: End of file

  Error termination. Backtrace:
  #0  0x7f4a90f4dedb in next_record_r
    at ../../../libgfortran/io/transfer.c:3144
  #1  0x7f4a90f50327 in finalize_transfer
    at ../../../libgfortran/io/transfer.c:3629
  #2  0x418e42 in ???
  #3  0x4ae5d0 in ???
  #4  0x495ace in ???
  #5  0x4060bc in ???
  #6  0x7f4a907bf554 in ???
  #7  0x40a2c6 in ???
  #8  0xffffffffffffffff in ???
  ```

<!-- ROADMAP -->

## Roadmap

See the [open issues](https://github.com/igordub/md-gst/issues) for a list of proposed features (and known issues).

<!-- CONTRIBUTING -->

## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- LICENSE -->

## License

Distributed under the MIT License. See `LICENSE` for more information.

<!-- CONTACT -->

## Contact

Igors Dubanevics - [igordub](https://github.com/igordub)

Project Link: [https://github.com/igordub/md-gst](https://github.com/igordub/md-cap)

<!-- ACKNOWLEDGEMENTS -->

## Acknowledgements

Thanks to [Dr Sarah Harris](https://astbury.leeds.ac.uk/people/dr-sarah-harris/) and [Dr Geoffery Wells](https://www.ucl.ac.uk/pharmacy/people/dr-geoffrey-wells) for their patient and encouraging support in MD simulation prepartaion and analysis.

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

<!-- APENDIX -->

## Apendix

Scritps for the AMBER simulations were taken from [David Burnell's PhD thesis](http://etheses.dur.ac.uk/11055) (Appendix A).

### A.1 Energy Minimisation

#### A.1.1 Minimising Solvent

```
CAP_2cAMP: initial minimisation solvent + ions
&cntrl
imin = 1, maxcyc = 10000,
ncyc = 5000, ntb = 1,
ntr = 1, cut = 10.0,
ntwx = 100
/
Hold the Protein fixed
500.0
RES 1 403
END
END
```

#### A.1.2 Minimising Solute

```
CAP_2cAMP: initial minimisation whole system
&cntrl
imin = 1, maxcyc = 50000,
ncyc = 25000, ntb = 1,
ntr = 0, cut = 10.0
/
```

### A.2 Equilibration

#### A.2.1 Temperature Equilibration

```
CAP_2cAMP: heat equilibration
&cntrl
imin=0, irest=0,
nstlim=100000, dt=0.002,
ntc=2, ntf=2,
cut=10.0, ntb=1,
ntpr=500, ntwx=5000,
ntt=3, gamma_ln=1.0,
ntx=1, ig=-1,
tempi=0.0, temp0=300.0,
ntr=1, ioutfm=1
/
Keep CAP fixed with weak restraints
2.5
RES 1 403
END
END
```

#### A.2.2 Pressure Equilibration

```
CAP_2cAMP: density equilibration
&cntrl
imin=0, irest=1,
nstlim=25000, dt=0.002,
ntc=2, ntf=2,
ntx=5, taup=1.0,
cut=8.0, ntb=2,
ntpr=500, ntwx=500,
ntt=3, gamma_ln=2.0,
temp0=300.0, ig=-1,
ntr=1, ioutfm=1,
ntp=1
/
Keep CAP fixed with weak restraints
10.0
RES 1 403
END
END
```

### A.3 Production MD

```
CAP_2cAMP: 4000ps of production MD
&cntrl
imin = 0, irest = 1,
ntb = 2, pres0 = 1.0,
taup = 2.0, iwrap=1,
cut = 10.0, ntr = 0,
ntc = 2, ntf = 2,
temp0 = 300.0, ntx = 5,
ntt = 3, gamma_ln = 1.0,
ntp = 1, ig=-1,
nstlim = 2000000, dt = 0.002,
ntpr = 5000, ntwx = 5000,
ntwr = 10000, ioutfm=1
/
```

### A.4 Sample Script for MD

```bash
module purge
module load dot
module load amber/cuda/SPDP/gcc/12.0
PROTEIN=cap
VAR=2CAMP
#Execute Commands
pmemd.cuda -O -i ../../wat_min1.in -o ${PROTEIN}_${VAR}_min1.out -p ${PROTEIN}_${VAR}.prmtop -c
${PROTEIN}_${VAR}.inpcrd -r ${PROTEIN}_${VAR}_min1.rst -ref ${PROTEIN}_${VAR}.inpcrd
#
pmemd.cuda -O -i ../../wat_min2.in -o ${PROTEIN}_${VAR}_min2.out -p ${PROTEIN}_${VAR}.prmtop -c
${PROTEIN}_${VAR}_min1.rst -r ${PROTEIN}_${VAR}_min2.rst
#
pmemd.cuda -O -i ../../heat_nc.in -o ${PROTEIN}_${VAR}_heat.out -p ${PROTEIN}_${VAR}.prmtop -c
${PROTEIN}_${VAR}_min2.rst -r ${PROTEIN}_${VAR}_heat.rst -ref ${PROTEIN}_${VAR}_min2.rst -x
${PROTEIN}_${VAR}_heat.nc
#
pmemd.cuda -O -i ../../density_nc.in -o ${PROTEIN}_${VAR}_density.out -p
${PROTEIN}_${VAR}.prmtop -c ${PROTEIN}_${VAR}_heat.rst -r ${PROTEIN}_${VAR}_density.rst -ref
${PROTEIN}_${VAR}_heat.rst -x ${PROTEIN}_${VAR}_density.nc
#
pmemd.cuda -O -i ../../md-prod_nc.in -o ${PROTEIN}_${VAR}_mdprod1.out -p
${PROTEIN}_${VAR}.prmtop -c ${PROTEIN}_${VAR}_density.rst -r ${PROTEIN}_${VAR}_mdprod1.rst -x
${PROTEIN}_${VAR}_mdprod1.nc
```

<!-- REFERENCES -->

## References

- [David Burnell's PhD Thesis](http://etheses.dur.ac.uk/11055)
- [Running Minimization and MD (in explicit solvent)](https://ambermd.org/tutorials/basic/tutorial1/section5.htm)
- [Simulating a pharmaceutical compound using antechamber and the Generalized Amber Force Field](https://ambermd.org/tutorials/basic/tutorial4b/index.php)

<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
-->

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
