#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 09:32:34 2023

@author: mojtaba
"""


from ccdc import io
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import rc
from ase.io import read
import spglib
# from tqdm import tqdm
# import molpack_enkelt as mp

path = '/home/mojtaba/Dropbox/PhD/dft_project/csp/ML_relax/systems/MOGBOB/minus_1'



# =============================================================================
# CSP part
# # =============================================================================
csp_DF = pd.read_csv(path + '/structures.csv')

cutoff = 2000

csp_energy = csp_DF['energy'].to_numpy()[:cutoff]
csp_energy -= min(csp_energy)
csp_density = csp_DF['density'].to_numpy()[:cutoff]
csp_spg = csp_DF['spacegroup'].to_numpy()[:cutoff]
csp_id = csp_DF['id'].to_list()[:cutoff]

# similar_structures = ['napht.dioxy-QR-1-9993-3', 'dioxy.napht-QR-1-19978-3' ]

# for STR in similar_structures:
    # if STR not in csp_id:
        # csp_id.append(STR)


for i in csp_id:
    crystal_reader = io.CrystalReader(path + '/structures/'+ i +'.cif')
    with io.CrystalWriter(path +'/ML/structures/'+i+'.cif', append=True) as crystal_writer:
            crystal_writer.write(crystal_reader[0])

# Dataset = []
# Dataset += csp_id[:11]

# LEN_atoms = []
# for system in tqdm(range(len(csp_id)), desc='Reading crystal structures', colour='red'):
    # atoms = read(path + '/structures/csd_corrected/' + csp_id[system] + '.cif')
    # LEN_atoms.append(len(atoms))

# max_len = np.max(LEN_atoms)
# csp_id_max_len = csp_id[np.array(LEN_atoms) == max_len]
    
    
    
