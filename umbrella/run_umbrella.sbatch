#!/bin/bash
# email on start, end, and abortion

#SBATCH --partition=andrewferguson-gpu
#SBATCH -o o.out
#SBATCH -e e1.out
#SBATCH --nodes=1            # SET NUM NODES
#SBATCH --ntasks-per-node=10  # SETS NUM MPI RANKS (1 PER GPU)
#SBATCH --cpus-per-task=1    # SET NUM THREADS
#SBATCH --account=pi-andrewferguson
#SBATCH --mem=8GB
#SBATCH --gres=gpu:1

folders=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09")
f_str="00 01 02 03 04 05 06 07 08 09"
NP=10

for folder in "${folders[@]}"; do
    cd "$folder"
    rm -f equil* umbrella*
    cp ../topol.top .
    gmx_mpi grompp -f ../npt_umbrella.mdp -c npt.gro -p topol.top -o equil.tpr -maxwarn 2
    cd ..
done
mpirun -np $NP gmx_mpi mdrun -deffnm equil -multidir $f_str -v -plumed plumed.dat -ntomp 1
for folder in "${folders[@]}"; do
    cd "$folder"
    gmx_mpi grompp -f ../npt_run.mdp -c equil.gro -t equil.cpt -p topol.top -o umbrella.tpr -maxwarn 2
    cd ..
done
mpirun -np $NP gmx_mpi mdrun -deffnm umbrella -multidir $f_str -v -plumed plumed.dat -ntomp 1

