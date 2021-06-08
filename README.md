# Neural Network Molecular Dynamics Kit

NNMDKit generates polymer topology from SMILES strings and automatically builds LAMMPS data and input files for molecular dynamics simulations

## Installation

```bash
pip install .
```

### Requirements
* [EMC](http://montecarlo.sourceforge.net/emc/Welcome.html)
* [RDKit](https://www.rdkit.org/)

Note that, both [EMC](http://montecarlo.sourceforge.net/emc/Welcome.html) and [RDKit](https://www.rdkit.org/) are required to be installed manually.
NNMDKit requires EMC to create polymer structures. To configure the integration with EMC, two environment variables are required to be addded to locate 
your EMC executable (`emc_linux64` for Linux, `emc_macos` for MacOS, or `emc_win32` for Windows) and setup tool (`emc_setup.pl`). Add the paths of the 
EMC executable and setup tool as environment variables "EMC_EXEC" and "EMC_SETUP", respectively.

For Ramprasad Group users on Tyrion2, simply type in the command line:

```
echo "export EMC_EXEC=/data/kevin/EMC/bin/emc_linux64" >> ~/.bashrc
echo "export EMC_SETUP=/data/kevin/EMC/scripts/emc_setup.pl" >> ~/.bashrc
```

## Example
Below is an example python script where we use NNMDKit to generate LAMMPS data and input files for Tg measurement with a list of SMILES strings.
```python
from nnmdkit.mkinput import write_data, write_lammps_input, write_pbs

smiles = ['*CC*', '*CC(*)C', '*CC(*)CC', '*CC(*)CCC', '*CC(*)CCCC','*CC(*)c1ccccc1']

for s in smiles:
    write_data(smiles=s, mw=10000, ntotal=3000, density=0.5, output_dir=s)

    write_lammps_input(Tinit=600, Tfinal=100, Tinterval=25, step=1000000, 
                       NN_POTENTIAL='~/p-rramprasad3-0/NNMD/potential_saved', output_dir=s)
                       
    write_pbs(jobname=s, nodes=2, ppn=24, walltime='48:00:00',
              LAMMPS_EXEC='~/p-rramprasad3-0/NNMD/lmp', output_dir=s)
```

## License & copyright
Ramprasad Group, Georgia Tech, USA

Licensed under the [MIT License](LICENSE). 