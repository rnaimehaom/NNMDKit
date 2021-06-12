from nnmdkit.util import Util

class Job:
    def __init__(self, jobname, project, nodes, ppn, walltime, LAMMPS_EXEC):
        self.jobname = jobname
        self.project = project
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.LAMMPS_EXEC = LAMMPS_EXEC

    def write_pbs(self, output_dir):
        Util.build_dir(output_dir)
        pbs_fname = 'job.pbs'
        with open(output_dir + '/' + pbs_fname, 'w') as f:
            f.write('#PBS -A {}\n'.format(self.project))
            f.write('#PBS -q inferno\n')
            f.write('#PBS -N {}\n'.format(self.jobname))
            f.write('#PBS -l nodes={}:ppn={}\n'.format(self.nodes, self.ppn))
            f.write('#PBS -l walltime={}\n'.format(self.walltime))
            f.write('#PBS -j oe\n')
            f.write('#PBS -o out.$PBS_JOBID\n')
            f.write('\n')
            f.write('cd $PBS_O_WORKDIR\n')
            f.write('mpirun -np {} {} -in lmp.in\n'.format(int(self.nodes*self.ppn), self.LAMMPS_EXEC))