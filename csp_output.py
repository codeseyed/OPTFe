#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 10:55:41 2022

@author: mojtaba
"""
# from ccdc import io
import pandas as pd
# import gemmi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from ase.io import read
# from matplotlib.ticker import ScalarFormatter  
import os
import spglib
import re
from molcrys import decompose
from ase.visualize import view
from matplotlib.ticker import FormatStrFormatter

csp_path = '/home/mojtaba/Dropbox/PhD/dft_project/csp/Design/'

# csp_path = '/home/mojtaba/Dropbox/PhD/dft_project/csp/ML_relax/systems/MAMPUM03/'


# # =============================================================================
# # correcting cif files with CSD
# # =============================================================================
# # for i in csp_id:
# #     crystal_reader = io.CrystalReader(csp_path + '/structures/'+ i +'.cif')
# #     with io.CrystalWriter(csp_path +'/structures/'+ i.replace(".", "") +'.cif', append=True) as crystal_writer:
# #             crystal_writer.write(crystal_reader[0])

# # csp_DF['id'] = [i.replace(".", "") for i in csp_id]


# # csp_to_molpack = csp_DF[['id', 'energy', 'density']].to_numpy()[:cutoff]
# # np.savetxt(csp_path + '/molpack_input.txt', csp_to_molpack, fmt='%s')

polar_spgs = [1, 3, 4, 5, 6, 7, 8, 9, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
            36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 75, 76, 77, 78, 79, 80,
            99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 143, 144,
            145, 146, 157, 156, 158, 159, 160, 161, 168, 169, 170, 171, 172, 173,
            183, 184, 185, 186]

# # sg = gemmi.SpaceGroup('61')
# # print(sg.short_name())

# # csp_spg_HM_notation = [gemmi.SpaceGroup(str(i)).short_name() for i in csp_spg]

# # polar_mask = [True if i in polar_spgs else False for i in csp_spg ]
# # nonpolar_mask = [False if i in polar_spgs else True for i in csp_spg ]


# # =============================================================================
# # DFT part
# # =============================================================================

combtion = 'e11'

csp_density_dft = []
csp_energy_dft = []
csp_spgs_dft = []
# # # forces = []
# # csp_removed_id = []
CSP_ID = []

os.chdir(csp_path + '/full_CSP/'+combtion+'/polar/structures_short/')
csp_id_1 = os.listdir('.')

# csp_DF_1 = pd.read_csv(csp_path + '/full_CSP/'+combtion+'/polar/structures.csv')
# cutoff = 1000

# csp_energy_dma = csp_DF_1['energy'].to_numpy()[:cutoff]
# # # csp_energy -= min(csp_energy)
# # # csp_density = csp_DF_1['density'].to_list()[:cutoff]
# # # csp_spg = csp_DF_1['spacegroup'].to_list()[:cutoff]
# # csp_id_1 = csp_DF_1['id'].to_list()[:cutoff]

for l,i in enumerate(csp_id_1):
    
    
    try:
        
        atoms = read(csp_path + '/full_CSP/'+combtion+'/polar/structures_short/'+ i + '/CONTCAR')
        outcar = read(csp_path + '/full_CSP/'+combtion+'/polar/structures_short/'+ i + '/OUTCAR')
        energy = outcar.get_total_energy()
        mols = decompose(atoms)
        
        # mol_1 = read(csp_path + '/MOLS/' + combtion[0] + '/minus_1/OUTCAR').get_total_energy()
        # mol_2 = read(csp_path + '/MOLS/' + combtion[1] + '/minus_1/OUTCAR').get_total_energy()
        
        # num_mols = len(mols)/2
        
        
        # force = outcar.get_forces()
        # forces.append(np.sum(np.linalg.norm(force, axis=1)))
        energy /= len(mols)
        # coh_energi = (energy - num_mols*mol_1 - num_mols*mol_2)/(len(mols))
    # energy in Kj/mol
        # energy *= 96.486
    # density in g/cm3
        density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
        csp_density_dft.append(density)
        csp_energy_dft.append(energy)
        space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
        csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
        CSP_ID.append(i)
    except Exception as e:
            print(e)
            print(i)
    
    if i == 'napht.dioxy-QR-1-10828-3':
        print(l)
                
        atoms = read(csp_path + '/full_CSP/'+combtion+'/polar/structures_short/'+ i + '/double_neg/CONTCAR')
        outcar = read(csp_path + '/full_CSP/'+combtion+'/polar/structures_short/'+ i + '/double_neg/OUTCAR')
        energy = outcar.get_total_energy()
        mols = decompose(atoms)
        
        # mol_1 = read(csp_path + '/MOLS/' + combtion[0] + '/minus_1/OUTCAR').get_total_energy()
        # mol_2 = read(csp_path + '/MOLS/' + combtion[1] + '/minus_1/OUTCAR').get_total_energy()
        
        # num_mols = len(mols)/2
        
        
        # force = outcar.get_forces()
        # forces.append(np.sum(np.linalg.norm(force, axis=1)))
        energy /= len(mols)
        # coh_energi = (energy - num_mols*mol_1 - num_mols*mol_2)/(len(mols))
    # energy in Kj/mol
        # energy *= 96.486
    # density in g/cm3
        density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
        csp_density_dft.append(density)
        csp_energy_dft.append(energy)
        space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
        csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
        CSP_ID.append(i)
        
        # csp_removed_id.append(i)
        # atoms = read(csp_path + '/minus_1/structures/corrected/'+ i + '/POSCAR')
    # density in g/cm3
        # density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
        # csp_density_dft.append(density)
        # csp_energy_dft.append(10.**6)
        # space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
        # csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
        # print(i)
        # csp_removed_id.append(i)


# csp_density_dft.append(1.34474376)
# csp_energy_dft.append(-5.87221836125)
# csp_spgs_dft.append(2)
# CSP_ID.append('napht.dioxy-QR-4-42001-3')

os.chdir(csp_path + '/full_CSP/'+combtion+'/para/structures_short/')
csp_id_2 = os.listdir('.')

# # # # csp_DF_2 = pd.read_csv(csp_path + 'minus_2/structures.csv')
# # # # # cutoff = 1000

# # # # # csp_id_2 = csp_DF_2['id'].to_list()[:cutoff]


for i in csp_id_2:
    
    
    try:
        
        atoms = read(csp_path + '/full_CSP/'+combtion+'/para/structures_short/'+ i + '/CONTCAR')
        outcar = read(csp_path + '/full_CSP/'+combtion+'/para/structures_short/'+ i + '/OUTCAR')
        energy = outcar.get_total_energy()
        mols = decompose(atoms)
        # force = outcar.get_forces()
        # forces.append(np.sum(np.linalg.norm(force, axis=1)))
        # energy /= len(atoms)
        # mol_1 = read(csp_path + '/MOLS/' + combtion[0] + '/neutral/OUTCAR').get_total_energy()
        # mol_2 = read(csp_path + '/MOLS/' + combtion[1] + '/neutral/OUTCAR').get_total_energy()
        
        # num_mols = len(mols)/2
        
        
        # force = outcar.get_forces()
        # forces.append(np.sum(np.linalg.norm(force, axis=1)))
        energy /= len(mols)
        # coh_energi = (energy - num_mols*mol_1 - num_mols*mol_2)/(len(mols))
    # energy in Kj/mol
        # energy *= 96.486
    # density in g/cm3
        density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
        csp_density_dft.append(density)
        csp_energy_dft.append(energy)
        space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
        csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
        CSP_ID.append(i)
    except Exception as e:
              print(e)
              print(i)
# #         # csp_id_2.remove(i)
# #         # atoms = read(csp_path + '/minus_2/structures/corrected/'+ i + '/POSCAR')
# #     # density in g/cm3
# #         # density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
# #         # csp_density_dft.append(density)
# #         # csp_energy_dft.append(10.**6)
# #         # space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
# #         # csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
        

# os.chdir(csp_path + '/neutral/structures_short/corrected')
# csp_id_3 = os.listdir('.')
# # # # csp_DF_3 = pd.read_csv(csp_path + 'neutral/structures.csv')
# # # # cutoff = 1000

# # # # csp_id_3 = csp_DF_3['id'].to_list()[:cutoff]


# for i in csp_id_3:
    
    
#     try:
        
#         atoms = read(csp_path + 'neutral/structures_short/corrected/'+ i + '/CONTCAR')
#         outcar = read(csp_path + 'neutral/structures_short/corrected/'+ i + '/OUTCAR')
#         energy = outcar.get_total_energy()
#         # force = outcar.get_forces()
#         # forces.append(np.sum(np.linalg.norm(force, axis=1)))
#         energy /= len(atoms)
#     # energy in Kj/mol
#         # energy *= 96.486
#     # density in g/cm3
#         density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
#         csp_density_dft.append(density)
#         csp_energy_dft.append(energy)
#         space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
#         csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
#         CSP_ID.append(i)
#     except Exception as e:
#         print(e)

#         print(i)
        # csp_id_3.remove(i)
        # atoms = read(csp_path + '/neutral/structures/corrected/'+ i + '/POSCAR')
    # density in g/cm3
        # density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
        # csp_density_dft.append(density)
        # csp_energy_dft.append(10.**6)
        # space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
        # csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))


    #     atoms = read(csp_path + 'structures/corrected_sp/'+ i + '/POSCAR')
    #     energy = 0.
    #     # force = outcar.get_forces()
    #     # forces.append(np.sum(np.linalg.norm(force, axis=1)))
    #     energy /= len(atoms)
    # # energy in Kj/mol
    #     energy *= 96.486
    # # density in g/cm3
    #     density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
    #     csp_density_dft.append(density)
    #     csp_energy_dft.append(energy)
    #     space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
    #     csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))




# csp_id = csp_id_1 + csp_id_2 + csp_id_3



# # csp_id = [i for i in csp_id if i not in csp_removed_id]


# =============================================================================
# THE template packing structures
# =============================================================================
# atoms = read(csp_path + '/template_packing/'+combtion+'/polar/CONTCAR')
# outcar = read(csp_path + '/template_packing/'+combtion+'/polar/OUTCAR')

# energy = outcar.get_total_energy()

# energy /= len(atoms)
# density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
# csp_density_dft.append(density)
# csp_energy_dft.append(energy)
# space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
# csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
# CSP_ID.append(i)

# atoms = read(csp_path + '/template_packing/'+combtion+'/centro/charged/CONTCAR')
# outcar = read(csp_path + '/template_packing/'+combtion+'/centro/charged/OUTCAR')

# energy = outcar.get_total_energy()

# energy /= len(atoms)
# density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
# csp_density_dft.append(density)
# csp_energy_dft.append(energy)
# space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
# csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
# CSP_ID.append(i)

# atoms = read(csp_path + '/template_packing/'+combtion+'/centro/neutral/CONTCAR')
# outcar = read(csp_path + '/template_packing/'+combtion+'/centro/neutral/OUTCAR')

# energy = outcar.get_total_energy()

# energy /= len(atoms)
# density = np.sum(atoms.get_masses())*1.66054/atoms.get_volume()
    
# csp_density_dft.append(density)
# csp_energy_dft.append(energy)
# space_g = str.split(spglib.get_spacegroup(atoms))[1]
        
# csp_spgs_dft.append(int(re.findall(r'-?\d+\.?\d*', space_g)[0]))
# CSP_ID.append(i)









csp_energy_dft = np.array(csp_energy_dft)
csp_energy_dft =  csp_energy_dft - min(csp_energy_dft)
# csp_energy_dft *= 30

csp_density_dft = np.array(csp_density_dft)

# 


polar_energy, polar_density, polar_spg, polar_ID =  zip(*sorted(zip(csp_energy_dft[:20], csp_density_dft[:20], csp_spgs_dft[:20],  CSP_ID[:20])))
para_energy, para_density, para_spg, para_ID =  zip(*sorted(zip(csp_energy_dft[20:], csp_density_dft[20:], csp_spgs_dft[20:],  CSP_ID[20:])))


polar_energy = np.array(polar_energy)
polar_density = np.array(polar_density)
para_energy = np.array(para_energy)
para_density = np.array(para_density)

a1_FE_pack = [True, True, True, False, True, True, True, True, True, False,
              True, True, True, True, True, True, True, True, True, False, True, 
              ] + [False]*20 
e1_FE_pack = [True]*16 + [False]*24 

e10_FE_pack = [True, True, False, False, True, False, False,
               True, False, False ] + [False]*30 
e11_FE_pack = [True, False, False, False, False, False, True, True, True,
               False ] + [False]*31
csp_energy_dft, csp_density_dft, CSP_ID =  zip(*sorted(zip(csp_energy_dft, csp_density_dft, CSP_ID)))
csp_density_dft = np.array(csp_density_dft)
csp_energy_dft = np.array(csp_energy_dft)

polar_mask = [True if i in polar_spgs else False for i in polar_spg ]
polar_mask_para = [True if i in polar_spgs else False for i in para_spg ]


# =============================================================================

rc('font',**{'family':'serif','serif':['Kerkis']})
rc('text', usetex=True)
plt.rcParams['figure.figsize'] = [6, 5]
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(dpi=400)

# ax.scatter(x=csp_density[polar_mask_dft],y=csp_energy[polar_mask_dft], c='C1',s=15,  label='Polar spacegroups')
# ax.scatter(x=csp_density[nonpolar_mask_dft],y=csp_energy[nonpolar_mask_dft], c='C4',s=15)

           
ax.scatter(x=polar_density,y=polar_energy,
            c='#0b5394ff',s=45, alpha=0.5, label='Monovalent salt')

ax.scatter(x=para_density,y=para_energy,
            c='#e69138ff',s=45, alpha=0.5, label='Neutral co-crystal')
           
ax.scatter(x=polar_density[polar_mask],y=polar_energy[polar_mask]
            ,s=45, facecolors='none', edgecolors='black', label='Polar spacegroup')
ax.scatter(x=para_density[polar_mask_para],y=para_energy[polar_mask_para]
            ,s=45, facecolors='none', edgecolors='black')

ax.scatter(x=csp_density_dft[e11_FE_pack],y=csp_energy_dft[e11_FE_pack],facecolors='none', edgecolors='darkred',
            s=55,
            marker='D',   label='Connected PT path')
plt.annotate( 'Ferroelectric packing', xy=(1.3,-0.003),
              xytext=(1.29,0.005) , va = "top", ha="left"  ) 
# plt.annotate( 'Ferroelectric packing', xy=(1.721,-0.003),
              # xytext=(1.6,0.01) , va = "top", ha="left" , arrowprops=dict( arrowstyle="->",shrinkA=0 ) )
# FE_double_ind = [i for i in range(len(CSP_ID)) if CSP_ID[i] == 'napht.dioxy-QR-4-42001-3' ]
ax.scatter(x=csp_density_dft[20],y=csp_energy_dft[20],facecolors='darkviolet',
            s=35,
                label='Divalent salt')


ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax.set_ylabel(r' Relative energy ($ \rm eV$)')
ax.set_xlabel(r'Density ($\rm g/ cm^{3}$)')

ymin, ymax = ax.get_ylim()
xmin, xmax = ax.get_xlim()

# ax.scatter(x=1.44,y=100., c='C4',s=30, 
            # edgecolors='black',   label='Anti-FE packing')
# ax.scatter(x=1.44,y=100., c='#3d85c6ff',s=30, 
            # edgecolors='red',   label='FE packing')

ax.legend(loc='best', ncol=1)
# plt.legend(loc="upper right", bbox_to_anchor=(0.53,0.85))
# plt.ylim(ymin, ymax)
# plt.xlim(1.625, xmax)      

plt.show()     



# rc('font',**{'family':'serif','serif':['Kerkis']})
# rc('text', usetex=True)

# plt.rcParams.update({'font.size': 12})
# fig, ax = plt.subplots(dpi=400)

           
# ax.scatter(x=range(1, len(csp_energi_dev)+1),y=csp_energi_dev,
#             c='C0',s=25, alpha=0.4, label='vdW-DF2')
# ax.plot(range(1, len(csp_energi_dev)+1),csp_energi_dev,
#             c='C0',lw=1.)

# ax.scatter(x=range(1, len(csp_energi_dev)+1),y=csp_energy_dma,
#             c='C1',s=25, alpha=0.7, label='DMACRYS')

# ax.plot(range(1, len(csp_energi_dev)+1),csp_energy_dma,
#                     c='C1',lw=1.)   

# ax.set_ylabel(r' Relative energy ($ \rm eV$)')
# ax.set_xlabel(r'Structures ranking')

# ymin, ymax = ax.get_ylim()
# xmin, xmax = ax.get_xlim()


# ax.legend(loc='best', ncol=1, fontsize=12)

# ax.set_xticks(np.arange(1, 21))
# plt.show()     





                
