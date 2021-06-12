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
    