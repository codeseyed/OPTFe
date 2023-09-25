#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 14:15:32 2023

@author: Mojtaba & Elin
"""
import numpy as np
from ase.io import read, write
from ase.io.vasp import write_vasp_xdatcar
from ase.build.tools import sort
from ase.db import connect
from dscribe.descriptors import SOAP
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.decomposition import KernelPCA 
from kneed import KneeLocator
from tqdm import tqdm

# =============================================================================
# =============================================================================
def atoms_size_distribution(csp_id, path):
    lents = [len(read(path + '/ML/structures/' + i + '.cif')) for i in csp_id]
    unik_lents = np.unique(lents)
    frequency = [lents.count(i) for i in unik_lents]
    return dict(zip(unik_lents, frequency))

# =============================================================================
# =============================================================================
def read_database(database_name, atomslist, total, lim):
    print('---------------------------------------------------------------------')
    print('---------------------------------------------------------------------')
    database  = connect(database_name) 
    print('Database:', database_name)
    print('From index:', total)
    for i in tqdm(range(len(database)),desc='Reading crystal structures', colour='red'):
        if i > lim:
            total += 1
            atoms = atoms = database.get_atoms(id=1+i)
            atoms = sort(atoms)
            atomslist.append(atoms)
    print('To index:', total)
    return(atomslist, total)

# =============================================================================
# =============================================================================
def make_soap_descriptors(atomslist, total, elements):
    # Uncomment lines to use cartesian mesh in stead of atomic centres for soap
    # with positions_cart in stead of positions
    print('---------------------------------------------------------------------')
    print('---------------------------------------------------------------------')
    print('Make SOAP descriptors')
    X = []
    #grid_num = 3

    for i in tqdm(range(total),desc='Making descriptors', colour='red'):
        atoms = atomslist[i]   

        soaper = SOAP(
        rcut=10,
        nmax=12,
        lmax=12,
        species=elements,
        sparse=False,
        weighting={"function": "poly", "r0": 12, "m": 2, "c": 1, "d": 1},
        periodic=True,
        average='inner'
        )
    
        inds =  [i.index for i in atoms if i.symbol == 'N' or 'Cl']
    
        positions = atoms.positions #NB: positions for soap descriptors, not atoms!
        #cell = atoms.get_cell() 
        #range_xy = np.linspace(0, 1, grid_num)
        #x, y, z = np.meshgrid(range_xy, range_xy,  range_xy )
        #positions_scaled = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T
        #positions_cart = cell.cartesian_positions(positions_scaled)
        soap = soaper.create(atoms, positions=positions)
        X.append(soap)

    return(X)

# =============================================================================
# =============================================================================
def PCA_descriptors(X_new, plot=True):
    pca = KernelPCA(kernel='linear')
    X_kpca = pca.fit_transform(X_new)
    explained_variance = np.var(X_kpca, axis=0)
    VAR = explained_variance / np.sum(explained_variance)
    if plot:
        print('---------------------------------------------------------------------')
        print('---------------------------------------------------------------------')
        print('Making PCA convergence plot')

        plt.rcParams.update({'font.size': 12})
        fig, ax = plt.subplots(dpi=400)
        plt.plot(range(len(VAR)), VAR.cumsum(), marker='o', color='C3')
        ax.set_ylabel(r' Cumulative explianed variance')
        ax.set_xlabel(r'Number of PCA components')
        plt.savefig('PCA_components.png')

    pca = KernelPCA(n_components=40, kernel='linear')  # for our test case. 40 PCA components 
                                                        # with linear kernel manage to cluster accurately
                                                        
                                                        
    score_pca = pca.fit_transform(X_new)
    return(score_pca)

# =============================================================================
# =============================================================================
def cluster_structures(score_pca, plot=True):
    kmeans_kwargs = {'init': 'random','n_init': 10, 'max_iter': 500,
                 'random_state': 42}

    sse = []
    for k in range(1, 31):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(score_pca)
        sse.append(kmeans.inertia_)
    if plot:
        print('---------------------------------------------------------------------')
        print('---------------------------------------------------------------------')
        print('Making cluster convergence plot')
    
        plt.rcParams.update({'font.size': 12})
        fig, ax = plt.subplots(dpi=400)
        plt.plot(range(len(sse)), sse, marker='o', color='C1')
        ax.set_ylabel(r'Sum of squared errors')
        ax.set_xlabel(r'Number of Kmean clusters')
        #plt.show()
        plt.savefig('clusters_opt.png')

    kl = KneeLocator(range(1, 31), sse, curve="convex", direction="decreasing" )

    opt_num_cluster = kl.knee

    print('---------------------------------------------------------------------')
    print('---------------------------------------------------------------------')
    print('The optimum number of clusters are: ', opt_num_cluster)
    # fitting k means from transformed coordinates checking for optimum number
    # of clusters
    # Split into different cluster with K-means
    model = KMeans(n_clusters=opt_num_cluster, **kmeans_kwargs)
    model.fit(score_pca)
    labels = model.labels_
    return(labels, opt_num_cluster)

# =============================================================================
# =============================================================================
def select_random_structures_from_cluster(num_struc_per_cluster,labels, atomslist, opt_num_cluster, total_struct, Write=True):
    Random_selected_strs = []
    total_list = np.linspace(0, total_struct-1, total_struct)
    csp_id = np.array(total_list)

    selected_atoms = []

    for group in range(opt_num_cluster):
        cstr = csp_id[labels == group]
        A = np.random.choice(cstr, num_struc_per_cluster)
        Random_selected_strs += A.tolist()
    
    print('Selected structures')
    print(Random_selected_strs)
    
    for i in Random_selected_strs:
        selected_atoms.append(atomslist[int(i)])
    if write:
        write_vasp_xdatcar('XDATCAR_selected', selected_atoms)
    return(Random_selected_strs, selected_atoms)
    
    
# =============================================================================
# =============================================================================
print('---------------------------------------------------------------------')
print('---------------------------------------------------------------------')
print('SELECTING MD DATA FOR NEURAL-IL TRAINING')

# =============================================================================
# =============================================================================

db_name1 = "MOGBOB_CSP.json"

atomslist = []
total_struct = 0
lim = -1
# =============================================================================
# Read all .json files
# =============================================================================
at_list, total_struct = read_database(db_name1, atomslist, total_struct, lim)
# at_list, total_struct = read_database(db_name2, at_list, total_struct, lim)
# at_list, total_struct = read_database(db_name3, at_list, total_struct, lim)
# at_list, total_struct = read_database(db_name4, at_list, total_struct, lim)
# at_list, total_struct = read_database(db_name5, at_list, total_struct, lim)

elements = ['Br', 'C', 'H', 'N', 'O']

# =============================================================================
# Make descriptors for all timesteps in atomslist
# =============================================================================

X = make_soap_descriptors(at_list, total_struct, elements)
X_new = np.zeros((len(X), len(X[0])), dtype=float)    

for ind in tqdm(range(len(X)), desc='Making X', colour='red'):
    X_new[ind] = X[ind]

# =============================================================================
# PCA for SOAP descriptors 
# =============================================================================

score_pca = PCA_descriptors(X_new)

# =============================================================================
# Clustering PCs
# =============================================================================

labels, opt_num_cluster = cluster_structures(score_pca, plot=True)

# =============================================================================
# Selecting N structures randomly from each cluster 
# =============================================================================

# Random_selected_strs, selected_atoms = select_random_structures_from_cluster(50,labels, atomslist, opt_num_cluster, total_struct, Write=True)

# =============================================================================
# Write ASE .json with selected atoms objects for traning NeuralIL
# =============================================================================
print('---------------------------------------------------------------------')
print('---------------------------------------------------------------------')
print('Writing ASE .json file with selected structures')

# new_database = 'Random_structures_chosen_from_clusters'

# with connect('{0}.json'.format(new_database)) as database:
    # for i, atom in enumerate(selected_atoms):
        # database.write(atom)

        
# np.savetxt('Random_struct_from_clusters.txt', Random_selected_strs, fmt='%s')

