import numpy as np
import subprocess
import sys
import os

sim_num = sys.argv[1]

def _osremove(f):
    if os.path.isfile(f):
        os.remove(f)
        
def mdrun(deffnm, plumed=False, plumed_file='plumed.dat', np=1, nsteps=100000, checkpoint=False, checkpoint_file="md.cpt"):
    """
    Python wrapper for gmx mdrun -deffnm

    Parameters
    ----------
    deffnm : str
         File names for md run
    mpirun : bool
        Is this a multi-node run or not gmx (False) vs gmx_mpi (Default: True)
        number of processes (np)
    plumed : bool
        Turns plumed on/off
    np : int
        Number of processes to run mpirun on (Default 1 for non-mpi run)
    """
   
    
    if plumed:
        commands = ["gmx_mpi", "mdrun", "-deffnm", deffnm, "-ntomp", str(np),  "-nsteps", str(nsteps)]
    else:
        commands = ["gmx_mpi", "mdrun", "-deffnm", deffnm, "-ntomp", str(np), "-v", "-nsteps", str(nsteps)]

    if plumed:
        commands.extend(["-plumed", plumed_file])
    if checkpoint:
        commands.extend(["-cpi", checkpoint_file])
    subprocess.run(commands)
    return

first_run = True
while True:
    if first_run:
        mdrun(f"md{sim_num}", plumed=True, plumed_file=f"plumed{sim_num}.dat", np=os.environ.get("SLURM_CPUS_PER_TASK", 8), nsteps=1000000, checkpoint=False)
        first_run = False
    else:
        
        mdrun(f"md{sim_num}", plumed=True, plumed_file=f"plumed{sim_num}.dat", np=os.environ.get("SLURM_CPUS_PER_TASK", 8), nsteps=1000000, checkpoint=True, checkpoint_file=f"md{sim_num}.cpt")


    colvar = np.loadtxt(f"{sim_num}_annealing_output.txt")
    cur_dists = colvar[-1, 4:]
    final_rg0 = colvar[-10:, 1]
    final_rg1 = colvar[-10:, 2]
    final_rg2 = colvar[-10:, 3]
    if np.max(cur_dists) > 2:
        break
    if np.any([np.all(final_rg0 < 1.25), np.all(final_rg1 < 1.25), np.all(final_rg2 < 1.25)]):
        break

    _osremove(f"md{sim_num}.gro")
