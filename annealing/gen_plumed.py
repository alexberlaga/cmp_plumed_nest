import mdtraj as md
import os
import numpy as np  
import sys

file_num = sys.argv[1]

def stringlist(items):
    strlist = [str(a) for a in items]
    return ",".join(strlist)

md_top = md.load("npt.gro")
rg_atoms = np.array([md_top.top.select('resid 0 to 19 and name CA'), md_top.top.select('resid 20 to 39 and name CA'), md_top.top.select('resid 40 to 59 and name CA')]) + 1
mid_atoms = rg_atoms[:, 9]
pairs = [[mid_atoms[0], mid_atoms[1]], [mid_atoms[0], mid_atoms[2]], [mid_atoms[1], mid_atoms[2]]]

with open("many_cv_template.dat", "r") as f:
    lines = f.readlines()

rg_num = 0
dist_num = 0
dist_list = ["01", "12", "02"]

with open(f"plumed{file_num}.dat", "w") as f:
    for line in lines:
        if line.startswith("rg"):
            l = f"rg{rg_num}: GYRATION TYPE=RADIUS ATOMS=" + stringlist(rg_atoms[rg_num]) + "\n"           
            rg_num += 1
        elif line.startswith("dist"):
            l = f"dist{dist_list[dist_num]}: DISTANCE ATOMS=" + stringlist(pairs[dist_num]) + "\n"
            dist_num += 1
        elif line.startswith("PRINT"):
            l = f"PRINT ARG=rg0,rg1,rg2,dist01,dist12,dist02 STRIDE=10000 FILE={file_num}_annealing_output.txt\n"
        else:
            l = line
        f.write(l)
        
