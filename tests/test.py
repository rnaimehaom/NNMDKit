from nnmdkit.mkinput import write_data, write_lammps_input, write_pbs

smiles = ['*CC*', '*CC(*)C', '*CC(*)CC', '*CC(*)CCC', '*CC(*)CCCC','*CC(*)c1ccccc1']

for s in smiles:
    write_data(smiles=s, mw=10000, ntotal=3000, density=0.5, output_dir=s)
    write_lammps_input(Tinit=600, Tfinal=100, Tinterval=25, step=1000000, NN_POTENTIAL='~/p-rramprasad3-0/NNMD/potential_saved', output_dir=s)
    write_pbs(jobname=s, nodes=2, ppn=24, walltime='48:00:00', LAMMPS_EXEC='~/p-rramprasad3-0/NNMD/lmp', output_dir=s)