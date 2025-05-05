import MDAnalysis as mda
import matplotlib.pyplot as plt
import time
import pickle
import sys
import mdtraj as md
import os
import plumed
from MDAnalysis.analysis import rms
import pandas as pd
import pyemma.coordinates as coor
import deeptime
from mlcolvar.core.stats import TICA
import numpy as np
from sklearn.decomposition import PCA
from utils import *
from tqdm import tqdm
from mlcolvar.cvs import DeepTICA
from tqdm import tqdm

def compute_all_distances(traj):
    idxs = np.arange(traj.top.n_atoms)
    grid = np.array(np.meshgrid(idxs, idxs)).T.reshape(-1, 2)
    pairs = grid[grid[:, 0] > grid[:, 1]]
    dists = md.compute_distances(traj, pairs)
    return dists


# ----------------------------------------------------------------------------------------
# Fitting TICA on trajectories


tica = TICA(var_cutoff=1, lagtime=20)

for folder in ['EKGEKG', 'DMGCIG', 'PPGPPG', 'EDGDDG']:
    for i in range(1, 7):
        traj = md.load(f'{folder}/3_anneal/md{i}_whole.xtc', top=f'{folder}/3_anneal/md_whole.gro')
        n_frames = len(np.loadtxt(f'{folder}/3_anneal/{i}_annealing_output.txt'))
        
        traj_ca = traj.atom_slice(traj.top.select('name CA'))[:int(n_frames * 1.5)]
        
        dists = compute_all_distances(traj_ca)
        tica.partial_fit(dists)

with open("tica_short.pkl", "wb+") as f:
    pickle.dump(tica, f)


# ----------------------------------------------------------------------------------------
# Using MLColvar with TICA as CV


with open("tica_short.pkl", "rb") as f:
    tica = pickle.load(f)

import torch
import torch.nn as nn
import mlcolvar 
# Assume `tica_model` is your trained Deeptime TICA object
# Extract the transformation parameters
means = torch.from_numpy(tica.model.mean_0).to(torch.float32)
components = torch.from_numpy(tica.model.singular_vectors_left[:, 0].reshape((1,-1))).to(torch.float32)

class TICA_transform(mlcolvar.core.transform.Transform):
    def __init__(self, in_features, out_ffeatures):
        super(TICA_transform, self).__init__(in_features, out_features)
        self.fc1 = nn.Linear(len(means), 1, bias=False)
        self.fc1.weight = torch.nn.Parameter(components)
        
    def forward(self, X):
        x_centered = X.to(torch.float32) - means
        return self.fc1(x_centered)

tica_transform = TICA_transform(len(means), 1)


torch.save(tica_transform, 'torch_tica.pt')
# print("TICA output (PyTorch):", tica_output)


def make_plumed_file(folder):
    top = md.load(f'{folder}/3_anneal/npt.gro')
    idxs = top.top.select('name CA')
    grid = np.array(np.meshgrid(idxs, idxs)).T.reshape(-1, 2)
    pairs = grid[grid[:, 0] > grid[:, 1]]
    traj = md.load(f'{folder}/3_anneal/md0_whole.xtc', top=f'{folder}/3_anneal/md_whole.gro')
    traj_water = md.load(f'{folder}/3_anneal/md0.xtc', top=f'{folder}/3_anneal/npt.gro')
    traj_ca = traj.atom_slice(traj.top.select('name CA'))
    ca_dists = compute_all_distances(traj_ca)
    cvs = tica.transform(ca_dists)[:, 0]
    for tic_i in tqdm(np.arange(-1.0, 3.6, 0.1)):
        stri = "".join(str(np.round(tic_i, 1)).split('.')).replace('-', 'n').replace('n00', '00')
        if not os.path.isdir(f'{folder}/4_umbrella/{stri}'):
            os.mkdir(f'{folder}/4_umbrella/{stri}')
        shutil.copyfile('torch_tica.pt', f'{folder}/4_umbrella/{stri}/torch_tica.pt')
        frame = np.argmin(np.abs(cvs - tic_i))
        traj_water[frame].save_gro(f'{folder}/4_umbrella/{stri}/npt.gro')
        distance_string = ','.join([f'd{i+1}' for i in range(len(pairs))])
        with open(f'{folder}/4_umbrella/{stri}/plumed.dat', 'w') as f:
            for i, pair in enumerate(pairs):
                f.write(f'd{i+1}: DISTANCE ATOMS={pair[0]+1},{pair[1]+1}\n')
            f.write(f'tic: PYTORCH_MODEL FILE=torch_tica.pt ARG={distance_string}\n')
            f.write(f'u: RESTRAINT ARG=tic.node-0 KAPPA=1000.0 AT={tic_i}\n') 
            f.write(f'PRINT ARG=tic.node-0,u.bias FILE=COLVAR_{stri} STRIDE=500\n')

# make_plumed_file('EKGEKG')
# make_plumed_file('PPGPPG')
# make_plumed_file('DMGCIG')
# make_plumed_file('EDGDDG')
# make_plumed_file('EKGDKG')
# make_plumed_file('EKGEXG')
# make_plumed_file('EKGKEG')
# make_plumed_file('EXGEKG')