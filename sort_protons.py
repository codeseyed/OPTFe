#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:56:47 2021

@author: mojtaba
"""

# =============================================================================
# The script for identifying and movig the protons to switch to the other
# bistable state 
# =============================================================================

import numpy as np
from ase.io import read, write
from ase.build import separate
from itertools import product
import os

class Proton:

    def __init__(self, ase_atoms):
        self.ase_atoms = ase_atoms
        self.atoms_positions = ase_atoms.positions
        self.unit_vec = np.array(ase_atoms.cell)
        self.hetro_index = [i.index for i in ase_atoms if i.symbol == 'N' or i.symbol == 'O']
        self.H_index = [i.index for i in ase_atoms if i.symbol == 'H']
        
    
    def is_atom_inside(self, atom_coordin, structure):
        cell = structure.cell
        scaled_coordin = cell.scaled_positions(atom_coordin)
        return all( 0.0<i<1 for i in scaled_coordin) 



    def count(array, num):
        count = 0
        repeat_indeces = []   
        for i,j in enumerate(array):
            if j[0] == num:
                count += 1
                repeat_indeces.append(i)
        return repeat_indeces, count


    def in_the_mol(self, ind_i, ind_j):
        new_atoms = self.ase_atoms[[atom.index for atom in self.ase_atoms if atom.symbol !='H']]
        separated = separate(new_atoms)
        for mol in separated:
            mol_pos = mol.positions
            frst_cond = np.array([(np.round(i, decimals=2) == np.round(self.atoms_positions[ind_i],
                                decimals=2)).all() for i in mol_pos]).any()
            sec_cond = np.array([(np.round(i, decimals=2) == np.round(self.atoms_positions[ind_j],
                                decimals=2)).all() for i in mol_pos]).any()
            if frst_cond and sec_cond:
                yield True
            else:
                yield False
        
    def atom_periodic(self, min_index, atom_j):
        a = self.unit_vec[0]
        b = self.unit_vec[1]
        c = self.unit_vec[2]
        pbc_possibilities = [a*0, a, -a, b, -b, c, -c, a-b, a+b, b-a, -a-b, 
                             a-c, a+c, c-a, -a-c, b-c, b+c, c-b, -b-c, 
                             a+b+c, a+b-c, a-b+c, -a+b+c, -a-b+c, -a+b-c, a+-b-c,
                             -a-b-c]
        atom_j_pos = np.zeros((len(pbc_possibilities), 3), dtype=float)
        for i in range(len(pbc_possibilities)):
            atom_j_pos[i, :] = atom_j + pbc_possibilities[i] 
        return atom_j_pos[min_index]

    def min_dist(self, atom_i, atom_j):
        a = self.unit_vec[0]
        b = self.unit_vec[1]
        c = self.unit_vec[2]
        pbc_possibilities = [a*0, a, -a, b, -b, c, -c, a-b, a+b, b-a, -a-b, 
                             a-c, a+c, c-a, -a-c, b-c, b+c, c-b, -b-c, 
                             a+b+c, a+b-c, a-b+c, -a+b+c, -a-b+c, -a+b-c, a+-b-c,
                             -a-b-c]
        distances_array = np.zeros(len(pbc_possibilities), dtype=float)
        for i in range(len(pbc_possibilities)):
            distances_array[i] = np.sqrt(np.sum((atom_i - (atom_j + pbc_possibilities[i]))**2))
        
        min_index = np.argmin(distances_array)
        atom_j = self.atom_periodic(min_index ,atom_j)
        return (distances_array[min_index], atom_j, min_index, pbc_possibilities)



    def h_hetro_indeces(self):
        h_hetro_index = []
        for i,j in product(self.hetro_index, self.H_index):
            if  self.min_dist(self.atoms_positions[i], self.atoms_positions[j])[0] < 1.1:
                h_hetro_index.append((i,j))
    # element_counts = [h_hetro_index[count(i) for i in h_hetro_index[:,0]]]
        return np.array(h_hetro_index)
    

    def hetro_sort(self, ind):
        hetro_dis = []
        new_hetro = []
        for i,j in product([ind], self.hetro_index):
            if i != j and not np.array(list(self.in_the_mol(ind, j))).any():
                distance = self.min_dist(self.atoms_positions[i], self.atoms_positions[j])[0]
                 # set a distance cutoff for hydrogen trnasfer 
                if distance < 3.6:
                    hetro_dis.append(distance)
                    new_hetro.append(j)
        a, b = zip(*sorted(zip(hetro_dis, new_hetro)))
        return b



    def dot_prodcuts(self):
        h_hetro_index = self.h_hetro_indeces()
        
        dot_product = np.zeros((len(h_hetro_index), len(self.hetro_index) -1), dtype=float)
        dot_sort = []
        for i in range(len(h_hetro_index)):
            hetro_sorted = self.hetro_sort(h_hetro_index[i][0])
            for j in range(len(hetro_sorted)):
                h_direction = self.min_dist(self.atoms_positions[h_hetro_index[i][0]],
                                            self.atoms_positions[h_hetro_index[i][1]])[1] - self.atoms_positions[h_hetro_index[i][0]]
                hetro_direction = self.min_dist(self.atoms_positions[h_hetro_index[i][0]],
                                                self.atoms_positions[hetro_sorted[j]])[1] - self.atoms_positions[h_hetro_index[i][0]]
                dot_product[i][j] = np.dot(h_direction, hetro_direction)
                if j != 0 and abs(dot_product[i][j] - dot_product[i][j-1]) < 1. :
                    dot_sort.append(hetro_sorted[j-1])
                    break
            else:
                dot_sort.append(hetro_sorted[np.argmax(dot_product[i,:])])
 
    
        hetro_couple = [(i,j) for i,j in zip(h_hetro_index[:,0], dot_sort)] 
    
        return np.array(hetro_couple)
        

class NEB:
    def __init__(self, ase_initial_atoms, ase_final_atoms, neb_directory):
        
        self.ase_initial_atoms = ase_initial_atoms
        self.ase_final_atoms = ase_final_atoms
        self.neb_directory = neb_directory
        
        
    
    def is_atom_inside(self, atom_coordin, structure):
        cell = structure.cell
        scaled_coordin = cell.scaled_positions(atom_coordin)
        return all( 0<i<1 for i in scaled_coordin) 
       
    
    def get_interpolated_images(self, images):
        pos_init = self.ase_initial_atoms.get_positions()
        unit_vec = np.array(self.ase_initial_atoms.cell)
        pos_fin = self.ase_final_atoms.get_positions()
        start_copy = self.ase_final_atoms.copy()
   
        distances_matrix = np.zeros((7, len(pos_init)), dtype=float)
        distances_matrix[0, :] = np.sqrt(np.sum((pos_fin - pos_init)**2, axis=1))
        distances_matrix[1, :] = np.sqrt(np.sum((pos_fin - (pos_init - unit_vec[0]))**2, axis=1))
        distances_matrix[2, :] = np.sqrt(np.sum((pos_fin - (pos_init + unit_vec[0]))**2, axis=1))
        distances_matrix[3, :] = np.sqrt(np.sum((pos_fin - (pos_init - unit_vec[1]))**2, axis=1))                             
        distances_matrix[4, :] = np.sqrt(np.sum((pos_fin - (pos_init + unit_vec[1]))**2, axis=1))               
        distances_matrix[5, :] = np.sqrt(np.sum((pos_fin - (pos_init - unit_vec[2]))**2, axis=1))
        distances_matrix[6, :] = np.sqrt(np.sum((pos_fin - (pos_init + unit_vec[2]))**2, axis=1))
        min_indices = np.argmin(distances_matrix, axis=0)
        intepolated_img = np.zeros((images, len(pos_init), len(pos_init[0])), dtype=float)
        for index in range(len(min_indices)):
            if min_indices[index] == 0:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*pos_init[index] + j/(images+1)*pos_fin[index]
                    intepolated_img[j-1][index] = pos 
        
  
            elif min_indices[index] == 1:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*(pos_init[index] - unit_vec[0]) + j/(images+1)*pos_fin[index]
                    if self.is_atom_inside(pos, self.ase_final_atoms):
                        intepolated_img[j-1][index] = pos
                    else:
                        intepolated_img[j-1][index] = pos + unit_vec[0]

            elif min_indices[index] == 2:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*(pos_init[index] + unit_vec[0]) + j/(images+1)*pos_fin[index]
                    if self.is_atom_inside(pos, self.ase_initial_atoms):
                        intepolated_img[j-1][index] = pos
                    else:
                        intepolated_img[j-1][index] = pos - unit_vec[0]
                    
            elif min_indices[index] == 3:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*(pos_init[index] - unit_vec[1]) + j/(images+1)*pos_fin[index]
                    if self.is_atom_inside(pos, self.ase_initial_atoms):
                        intepolated_img[j-1][index] = pos
                    else:
                        intepolated_img[j-1][index] = pos + unit_vec[1]        

            elif min_indices[index] == 4:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*(pos_init[index] + unit_vec[1]) + j/(images+1)*pos_fin[index]
                    if self.is_atom_inside(pos, self.ase_initial_atoms):
                        intepolated_img[j-1][index] = pos
                    else:
                        intepolated_img[j-1][index] = pos - unit_vec[1] 
        
            elif min_indices[index] == 5:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*(pos_init[index] - unit_vec[2]) + j/(images+1)*pos_fin[index]
                    if self.is_atom_inside(pos, self.ase_initial_atoms):
                        intepolated_img[j-1][index] = pos
                    else:
                        intepolated_img[j-1][index] = pos + unit_vec[2]        

            elif min_indices[index] == 6:
                for j in range(1,images+1):
                    pos = (1- (j/(images+1)))*(pos_init[index] + unit_vec[2]) + j/(images+1)*pos_fin[index]
                    if self.is_atom_inside(pos, self.ase_initial_atoms):
                        intepolated_img[j-1][index] = pos
                    else:
                        intepolated_img[j-1][index] = pos - unit_vec[2] 
            
        for stru in range(1, images +1):
            start_copy.set_positions(intepolated_img[stru-1])
            if not os.path.exists(self.neb_directory + str(0)+str(stru)):
                os.mkdir(self.neb_directory + str(0)+str(stru))
            
            write(self.neb_directory+str(0)+str(stru)+'/POSCAR',
                        start_copy, format='vasp')
    def get_iamges_energies(self, images):
        energies = np.zeros(images+1, dtype=float)
        
        outcar = read(self.neb_directory + '00/OUTCAR')
        energies[0] = outcar.get_total_energy() 
        for i in range(1, images+1):
            outcar = read(self.neb_directory + '0' + str(i) + '/OUTCAR')
            energies[i] = outcar.get_total_energy() 
        # outcar = read(self.neb_directory + 'fin/OUTCAR')
        # energies[-1] = outcar.get_total_energy()
        return energies
        
        
       
                


            
    
            
        

        


            
        
        



