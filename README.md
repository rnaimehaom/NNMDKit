# Neural Network Molecular Dynamics Kit

NNMDKit generates polymer topology from SMILES strings and automatically builds LAMMPS data and input files for molecular dynamics simulations

## Installation

```
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