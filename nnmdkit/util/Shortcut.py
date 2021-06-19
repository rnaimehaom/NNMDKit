import os
import glob
from rdkit.Chem import Descriptors, MolFromSmiles
from subprocess import call


def write_data(smiles,
               mw,
               ntotal,
               density,
               output_dir,
               output_prefix='system',
               tmp_ff='opls-aa',
               terminator='*[H]',
               cleanup=True):

    try:
        os.mkdir(output_dir)
    except OSError:
        pass

    previous_dir = os.getcwd()
    os.chdir(output_dir)

    # Write .esh file required to run EMC
    tmp_eshfile = '{}.esh'.format(output_prefix)
    RU_mw = Descriptors.ExactMolWt(MolFromSmiles(smiles))
    chainlength = int(mw / RU_mw)
    with open(tmp_eshfile, 'w') as f:
        f.write('#!/usr/bin/env emc_setup.pl\n')
        f.write('ITEM OPTIONS\n')
        f.write('replace true\n')
        f.write('field {}\n'.format(tmp_ff))
        f.write('density {}\n'.format(density))
        f.write('ntotal {}\n'.format(ntotal))
        # f.write('emc_execute true\n')
        f.write('ITEM END\n')
        f.write('\n')
        f.write('ITEM GROUPS\n')
        f.write('RU {},1,RU:2\n'.format(smiles))
        f.write('terminator {},1,RU:1,1,RU:2\n'.format(terminator))
        f.write('ITEM END\n')
        f.write('\n')
        f.write('ITEM CLUSTERS\n')
        f.write('poly alternate 1\n')
        f.write('ITEM END\n')
        f.write('\n')
        f.write('ITEM POLYMERS\n')
        f.write('poly\n')
        f.write('1 RU,{},terminator,2\n'.format(chainlength))
        f.write('ITEM END\n')

    # Run EMC
    EMC_SETUP = os.environ.get('EMC_SETUP')
    EMC_EXEC = os.environ.get('EMC_EXEC')

    try:
        call('{} {}'.format(EMC_SETUP, tmp_eshfile), shell=True)
    except BaseException:
        print('problem setting up emc.')

    try:
        call('{} build.emc'.format(EMC_EXEC), shell=True)
    except BaseException:
        print('problem running emc.')

    # Read in the EMC generated data file for modification
    headerlines = []
    masslines = []
    atomlines = []
    readheader = True
    readmass = False
    readatom = False
    with open('{}.data'.format(output_prefix), 'rt') as lines:
        for line in lines:
            # Remove everything after the bond section in the data file
            if 'Bonds' in line:
                break
            if readheader and not any(
                    x in line
                    for x in ['bond', 'angle', 'dihedral', 'improper']):
                headerlines.append(line)
            if 'Atoms' in line:
                readmass = False
                readatom = True
            elif readmass and line != '\n':
                masslines.append(line)
            elif readatom and line != '\n':
                atomlines.append(line)
            elif 'Masses' in line:
                readheader = False
                readmass = True

    # Make each element the same type of particle (e.g. no difference between aromatic C and regular C)
    # Convert repeated type numbers to a unique type number
    unique_mass = []
    typeconvertdict = {}
    for line in masslines:
        fields = line.split()
        if fields[1] not in unique_mass:
            unique_mass.append(fields[1])
        typeconvertdict[fields[0]] = len(unique_mass)

    # Rewrite data output files with modified atom types and charges
    with open('{}.data'.format(output_prefix), 'wt') as out:
        for line in headerlines:
            if 'atom types' in line:
                out.write('{:>12}  atom types\n'.format(len(unique_mass)))
            else:
                out.write(line)
        out.write('\n')
        for n, i in enumerate(unique_mass):
            out.write('{:>8} {:>10}\n'.format(n + 1, i))
        out.write('\n')
        out.write('Atoms\n')
        out.write('\n')
        for line in atomlines:
            fields = line.split()
            new_atomtype = typeconvertdict[fields[2]]
            out.write(
                '{0:>8} {1:>7} {atomtype:>3} {charge:>7} {4:>14} {5:>14} {6:>14}\n'
                .format(*fields, atomtype=new_atomtype, charge=0))

    # Clean up all EMC generated files except for the data file
    if cleanup:
        fnames = ['build.emc']
        fnames += glob.glob('*.esh')
        fnames += glob.glob('*.gz')
        fnames += glob.glob('*.in')
        fnames += glob.glob('*.vmd')
        fnames += glob.glob('*.params')
        for fname in fnames:
            try:
                os.remove(fname)
            except BaseException:
                print('problem removing {} during cleanup'.format(fname))

    os.chdir(previous_dir)


def write_lammps_input(Tinit,
                       Tfinal,
                       Tinterval,
                       step,
                       NN_POTENTIAL,
                       output_dir,
                       datafile_prefix='system',
                       pair_style='nn',
                       element='C H'):

    # Setting up settings parameters
    neighbor_skin = 2.0
    neighbor_every = 1
    timestep = 1
    thermo = 1000

    settings_fname = 'system.in.settings'
    with open(output_dir + '/' + settings_fname, 'w') as f:
        f.write('{:<15} {}\n'.format('pair_style', pair_style))
        f.write('{:<15} * * {} {}\n'.format('pair_coeff', NN_POTENTIAL,
                                            element))
        f.write('\n')
        f.write('{:<15} {} bin\n'.format('neighbor', neighbor_skin))
        f.write('{:<15} delay 0 every {} check yes\n'.format(
            'neigh_modify', neighbor_every))
        f.write('\n')
        f.write(
            '{:<15} custom step temp density vol press ke pe ebond evdwl ecoul elong\n'
            .format('thermo_style'))
        f.write('{:<15} {}\n'.format('thermo', thermo))
        f.write('{:<15} {}\n'.format('timestep', timestep))

    # Setting up simulation parameters
    Pinit = 1
    Tmax = 1000  # 1000K
    Pmax = 49346.163  # 50000bar
    Tdamp = 100
    Pdamp = 1000

    # Setting up equilibration procedure: nvt/npt, length (ns), temp (K), pressure (bar)
    eq_step = [
        ['nvt', 50000, Tmax],
        ['nvt', 50000, Tinit],
        ['npt', 50000, Tinit, 0.02 * Pmax],
        ['nvt', 50000, Tmax],
        ['nvt', 100000, Tinit],
        ['npt', 50000, Tinit, 0.6 * Pmax],
        ['nvt', 50000, Tmax],
        ['nvt', 100000, Tinit],
        ['npt', 50000, Tinit, Pmax],
        ['nvt', 50000, Tmax],
        ['nvt', 100000, Tinit],
        ['npt', 5000, Tinit, 0.5 * Pmax],
        ['nvt', 5000, Tmax],
        ['nvt', 10000, Tinit],
        ['npt', 5000, Tinit, 0.1 * Pmax],
        ['nvt', 5000, Tmax],
        ['nvt', 10000, Tinit],
        ['npt', 5000, Tinit, 0.01 * Pmax],
        ['nvt', 5000, Tmax],
        ['nvt', 10000, Tinit],
        ['npt', 800000, Tinit, Pinit],
    ]
    totaltime_eq = 0
    for i in eq_step:
        totaltime_eq += i[1]

    # Writing LAMMPS input file
    lmp_input_fname = 'lmp.in'
    with open(output_dir + '/' + lmp_input_fname, 'w') as f:

        f.write(
            '# LAMMPS input file generated by high-throughput Polymer MD generator\n'
        )
        f.write('\n')

        f.write('### Initialization\n')
        f.write('{:<15} full\n'.format('atom_style'))
        f.write('{:<15} real\n'.format('units'))
        f.write('{:<15} {}.data\n'.format('read_data', datafile_prefix))
        f.write('{:<15} {}\n'.format('include', settings_fname))
        f.write('\n')
        f.write('\n')

        f.write('### Minimization\n')
        f.write('{:<15} cg\n'.format('min_style'))
        f.write('{:<15} 1.0e-6 1.0e-8 100000 10000000\n'.format('minimize'))
        f.write('{:<15} 0\n'.format('reset_timestep'))
        f.write('\n')
        f.write('\n')

        f.write('### Equilibration\n')
        f.write(
            '{:<15} dump1 all custom 10000 equil.lammpstrj id mol type q xs ys zs ix iy iz\n'
            .format('dump'))
        f.write('{:<15} {} equilibrated.restart\n'.format(
            'restart', totaltime_eq))
        f.write('\n')

        for n, i in enumerate(eq_step):
            if i[0] == 'nvt':
                f.write('{:<15} step{} all nvt temp {} {} {}\n'.format(
                    'fix', n + 1, i[2], i[2], Tdamp))
            elif i[0] == 'npt':
                f.write('{:<15} step{} all npt temp {} {} {} iso {} {} {}\n'.
                        format('fix', n + 1, i[2], i[2], Tdamp, i[3], i[3],
                               Pdamp))
            f.write('{:<15} {}\n'.format('run', i[1]))
            f.write('{:<15} step{}\n'.format('unfix', n + 1))
            f.write('\n')
        f.write('{:<15} dump1\n'.format('undump'))
        f.write('{:<15} 0\n'.format('reset_timestep'))
        f.write('\n')
        f.write('\n')

        f.write('### Production\n')
        f.write(
            '{:<15} dump2 all custom 10000 production.lammpstrj id mol type q xs ys zs ix iy iz\n'
            .format('dump'))
        f.write('{:<15} {} production.restart\n'.format('restart', step))
        f.write('{:<15} Rho equal density\n'.format('variable'))
        f.write('{:<15} Temp equal temp\n'.format('variable'))
        f.write(
            '{:<15} fDENS all ave/time {} {} {} v_Temp v_Rho file temp_vs_density\n'
            .format('fix', int(step / 100 / 4), 100, step))
        f.write('\n')

        f.write('{:<15} loop\n'.format('label'))
        f.write('{:<15} a loop {}\n'.format(
            'variable', int((Tinit - Tfinal) / Tinterval + 1)))
        f.write('{:<15} b equal {}-{}*($a-1)\n'.format('variable', Tinit,
                                                       Tinterval))
        f.write('{:<15} fNPT all npt temp $b $b {} iso {} {} {}\n'.format(
            'fix', Tdamp, Pinit, Pinit, Pdamp))
        f.write('{:<15} {}\n'.format('run', step))
        f.write('{:<15} fNPT\n'.format('unfix'))
        f.write('{:<15} a\n'.format('next'))
        f.write('{:<15} SELF loop\n'.format('jump'))
        f.write('{:<15} a delete\n'.format('variable'))


def write_pbs(jobname,
              nodes,
              ppn,
              walltime,
              LAMMPS_EXEC,
              output_dir,
              project='GT-rramprasad3-CODA20'):
    pbs_fname = 'job.pbs'
    with open(output_dir + '/' + pbs_fname, 'w') as f:
        f.write('#PBS -A {}\n'.format(project))
        f.write('#PBS -q inferno\n')
        f.write('#PBS -N {}\n'.format(jobname))
        f.write('#PBS -l nodes={}:ppn={}\n'.format(nodes, ppn))
        f.write('#PBS -l walltime={}\n'.format(walltime))
        f.write('#PBS -j oe\n')
        f.write('#PBS -o out.$PBS_JOBID\n')
        f.write('\n')
        f.write('cd $PBS_O_WORKDIR\n')
        f.write('mpirun -np {} {} -in lmp.in\n'.format(int(nodes * ppn),
                                                       LAMMPS_EXEC))


if __name__ == '__main__':
    smiles = [
        '*CC*', '*CC(*)C', '*CC(*)CC', '*CC(*)CCC', '*CC(*)CCCC',
        '*CC(*)c1ccccc1'
    ]

    for s in smiles:
        write_data(smiles=s, mw=10000, ntotal=3000, density=0.5, output_dir=s)
        write_lammps_input(
            Tinit=600,
            Tfinal=100,
            Tinterval=25,
            step=1000000,
            NN_POTENTIAL='~/p-rramprasad3-0/NNLMP/potential_saved',
            output_dir=s)
        write_pbs(jobname=s,
                  nodes=2,
                  ppn=24,
                  walltime='48:00:00',
                  LAMMPS_EXEC='~/p-rramprasad3-0/NNLMP/lmp',
                  output_dir=s)
