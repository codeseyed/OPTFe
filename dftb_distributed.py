#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:05:25 2022

@author: mojtaba
"""

from mpi4py import MPI
import numpy as np
import os
from shutil import copyfile
import pandas as pd

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

path = os.getcwd()

csp_DF = pd.read_csv(path + '/structures.csv')

cutoff = 20
csp_id = csp_DF['id'].to_list()[:cutoff]
csp_id.append('dioxy.napht-QR-1-19978-3')

if size < len(csp_id):
    
    step = len(csp_id)//size

    inds = np.arange(0, len(csp_id), step)

    if rank == (size -1):
        rank_systems = csp_id[inds[rank]:]
        for i in rank_systems:
            #os.mkdir(path + '/' + i)
            os.chdir(path + '/' + i)
            #copyfile(path + '/' + i + '.cif', path + '/' + i + '/' + i +'.cif')
            os.system('cspy-opt '+i+'.cif --dftb-calc --lattice-opt')

    else:
        rank_systems = csp_id[inds[rank]:inds[rank+1]]
        for i in rank_systems:
            #os.mkdir(path + '/' + i)
            os.chdir(path + '/' + i)
            #copyfile(path + '/' + i + '.cif', path + '/' + i + '/' + i +'.cif')
            os.system('cspy-opt '+i+'.cif --dftb-calc --lattice-opt')
            
else:
    
    for i in range(len(csp_id)):
        if rank == i:
            #os.mkdir(path + '/' + csp_id[i])
            os.chdir(path + '/' + csp_id[i])
            #copyfile(path + '/' + csp_id[i] + '.cif', path + '/' + csp_id[i] + '/' + csp_id[i] +'.cif')
            os.system('cspy-opt '+csp_id[i]+'.cif --dftb-calc --lattice-opt')
            
    
