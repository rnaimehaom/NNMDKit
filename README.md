# Neural Network Molecular Dynamics Kit

NNMDKit generates polymer topology from SMILES strings and automatically builds LAMMPS data and input files for molecular dynamics simulations

## Example
Below is an example python script where we use NNMDKit to generate LAMMPS data and input files for Tg measurement with a list of SMILES strings.
```python
import nnmdkit

smiles = ['*CC*', '*CC(*)C', '*CC(*)CC', '*CC(*)CCC', '*CC(*)CCCC','*CC(*)c1ccccc1']

for s in smiles:
    sys = nnmdkit.System(smiles=s, mw=10000, ntotal=3000, density=0.5)
    data = sys.write_data(output_dir=s)

    lmp = nnmdkit.Lammps(data, NN_POTENTIAL='~/p-rramprasad3-0/NNLMP/potential_saved')
    lmp.add_procedure('minimization', min_style='cg')
    lmp.add_procedure('equilibration', Tfinal=600, Pfinal=1, Tmax=800, Pmax=49346.163)
    lmp.add_procedure('Tg_measurement', Tinit=600, Tfinal=100, Tinterval=25, step=1000000)
    lmp.write_input(output_dir=s)
                       
    job = nnmdkit.Job(jobname=s, project='GT-rramprasad3-CODA20', nodes=2, ppn=24, 
                      walltime='48:00:00', LAMMPS_EXEC='~/p-rramprasad3-0/NNLMP/lmp')
    job.write_pbs(output_dir=s)
```
<p align="center">
    <img src='./tutorial/img/temp_vs_density.png' width="70%"/>
</p>

A tutorial on using NNMDKit to create simulations of hydrocarbon polymers can be found [here](https://github.com/Ramprasad-Group/NNMDKit/tree/master/tutorial/CaseStudies.Hydrobarbons.ipynb)

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

## License & copyright
Ramprasad Group, Georgia Tech, USA

Licensed under the [MIT License](LICENSE). 