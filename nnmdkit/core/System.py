import os
import glob
from nnmdkit.util import Util
from rdkit.Chem import Descriptors, MolFromSmiles
from subprocess import call


class System:
    '''nnmdkit.core.System.System

    Template object to contain class initialization settings

    Attributes:
        smiles: str
            SMILES string of the polymer (use * as connecting point)

        mw: int
            Molecular weight of the polymer

        ntotal: int
            Total number of atoms of the system

        density: float
            Density of the system
    '''
    def __init__(self, smiles, mw, ntotal, density):
        self.smiles = smiles
        self.mw = mw
        self.ntotal = ntotal
        self.density = density

    def write_data(self,
                   output_dir,
                   output_prefix='system',
                   tmp_ff='opls-aa',
                   terminator='*[H]',
                   cleanup=True):

        Util.build_dir(output_dir)

        previous_dir = os.getcwd()
        os.chdir(output_dir)

        # Write .esh file required to run EMC
        tmp_eshfile = '{}.esh'.format(output_prefix)
        RU_mw = Descriptors.ExactMolWt(MolFromSmiles(self.smiles))
        chainlength = int(self.mw / RU_mw)
        with open(tmp_eshfile, 'w') as f:
            f.write('#!/usr/bin/env emc_setup.pl\n')
            f.write('ITEM OPTIONS\n')
            f.write('replace true\n')
            f.write('field {}\n'.format(tmp_ff))
            f.write('density {}\n'.format(self.density))
            f.write('ntotal {}\n'.format(self.ntotal))
            # f.write('emc_execute true\n')
            f.write('ITEM END\n')
            f.write('\n')
            f.write('ITEM GROUPS\n')
            f.write('RU {},1,RU:2\n'.format(self.smiles))
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
        data_fname = '{}.data'.format(output_prefix)
        with open(data_fname, 'rt') as lines:
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
        with open(data_fname, 'wt') as out:
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
        return data_fname
