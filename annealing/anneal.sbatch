#!/bin/bash

# email on start, end, and abortion
#SBATCH -J EKGEKG.an.an$1

#SBATCH -o out.out
#SBATCH -e error.out
#SBATCH --partition=gm4-pmext
#SBATCH --qos=gm4
#SBATCH --nodes=1            # SET NUM NODES 
#SBATCH --ntasks-per-node=1  # SETS NUM MPI RANKS (1 PER GPU)
#SBATCH --cpus-per-task=10    # SET NUM THREADS 
#SBATCH --account=pi-andrewferguson
#SBATCH --mem=8GB
#SBATCH --gres=gpu:1

# THIS EXAMPLE USES 1 GPU NODE - 1 MPI TASK - 4 THREADS PER TASK

# SET NUMBER OF MPI TASKS 
# SET NUMBER OF MD STEPS

#LOAD GROMACS MODULE 
source /scratch/midway3/berlaga/miniconda3/bin/activate sams

NUMRUNS=$1
SIMNUM=$2

module load cuda/11.2

for i in $(seq $((SIMNUM)) $((SIMNUM+NUMRUNS-1)))
do
	echo $i
        python gen_plumed.py $i

        gmx_mpi grompp -f npt_run.mdp -c npt.gro -t npt.cpt  -p topol.top -o md${i}.tpr -maxwarn 2

        python run_sim.py $i
done

for i in $(seq $((SIMNUM)) $((SIMNUM+NUMRUNS-1)))
do
        echo 1 | gmx_mpi trjconv -f md${i}.xtc -o md${i}_whole.xtc -s md${i}.tpr -pbc whole

done
echo 1 | gmx_mpi trjconv -f md0.gro -o md_whole.gro -s md0.tpr -pbc whole
