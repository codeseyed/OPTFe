#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:59:43 2021

@author: mojtaba
"""

import numpy as np
import matplotlib.pyplot as plt
import GGAxclib as gga
import pandas as pd
import matplotlib as mpl
from matplotlib import rc

def pbe_corr_kernel(t,A):
  phi=1 #Spin neutral (for now)
  beta=0.066725
  gamma=0.031091
  # A = create_A(beta,gamma,phi)
  # A = 0.05
  At2 = A*t**2
  H = (gamma*phi**3)*np.log( 1 + beta/gamma*t**2*gga.aux_PBE(At2) )
  return H

def A(n):
    phi=1 #Spin neutral (for now)
    beta=0.066725
    gamma=0.031091

    return (beta/gamma)/(np.exp(-gga.eps_c_PW92_LDA_sn(n)/(phi**3*gamma)) -1)



def t(s,n):
    return 0.5*s*np.sqrt(gga.kf(n)*np.pi)

def categorizer(t_arr, n_arr, grid):
    n_space = np.linspace(0, 0.05, grid+1)
    t_space = np.linspace(0, 6, grid+1)
    t_n_matrix = np.zeros((grid, grid), dtype=float)
    count = 0
    for i in range(len(t_space)-1):
        count += 1
        print(count)
        for j in range(len(n_space)-1):
            
            my_index = np.where(( t_arr > t_space[i] ) & (t_arr < t_space[i+1]) \
                                                   & ( n_arr > n_space[j] ) & (n_arr < n_space[j+1]))
                
            t_n_matrix[i,j] = sum(n_arr[my_index[0]])
            
    return t_n_matrix

n_space = np.linspace(0, 0.05, 200+1)
t_space = np.linspace(0.06, 6, 200)

den_ts = abs(np.genfromtxt('/home/mojtaba/Desktop/dft_project/pt_article/systems/e/trans/pbe/e_ts_den.cube',
                        skip_header=16).flatten())
rdg_ts = np.genfromtxt('/home/mojtaba/Desktop/dft_project/pt_article/systems/e/trans/pbe/e_ts_rdg.cube',
                    skip_header=16).flatten()

den_min = abs(np.genfromtxt('/home/mojtaba/Desktop/dft_project/pt_article/systems/e/min/pbe/e_den_min.cube',
                        skip_header=16).flatten())
rdg_min = np.genfromtxt('/home/mojtaba/Desktop/dft_project/pt_article/systems/e/min/pbe/e_rdg_min.cube',
                    skip_header=16).flatten()


den_min = np.where((den_min <0.000001), 0.0, den_min)
rdg_min = np.where(den_min == 0.0, 0.0, rdg_min)


den_ts = np.where((den_ts <0.000001) , 0.0, den_ts)
rdg_ts = np.where(den_ts == 0.0, 0.0, rdg_ts)


trun_den_min = den_min[np.nonzero(den_min)[0]]
trun_den_ts = den_ts[np.nonzero(den_ts)[0]]
s_min = rdg_min[np.nonzero(rdg_min)[0]]
s_ts = rdg_ts[np.nonzero(rdg_ts)[0]]

# t_min = t(trun_rdg_min, trun_den_min)
# t_ts = t(trun_rdg_ts, trun_den_ts)


s_min = np.where((s_min > 6), 0.0, s_min)
s_ts = np.where((s_ts > 6), 0.0, s_ts)

trun_den_min = np.where(s_min == 0.0, 0.0, trun_den_min)
trun_den_ts = np.where(s_ts == 0.0, 0.0, trun_den_ts)

trun_den_min = trun_den_min[np.nonzero(trun_den_min)[0]]
trun_den_ts = trun_den_ts[np.nonzero(trun_den_ts)[0]]
s_min = s_min[np.nonzero(s_min)[0]]
s_ts = s_ts[np.nonzero(s_ts)[0]]


# t_n_arr_min = categorizer(s_min, trun_den_min, 200)
# t_n_arr_ts = categorizer(s_ts, trun_den_ts, 200)







# t_n_arr_min = np.genfromtxt('/home/mojtaba/Desktop/dft_project/pt_article/PBE_corr/t_n_array_min')
# t_n_arr_ts = np.genfromtxt('/home/mojtaba/Desktop/dft_project/pt_article/PBE_corr/t_n_array_ts')

# n = np.linspace(0, 1, 1000 )
# eps_pw92 = gga.eps_c_PW92_LDA_sn(n)

# eps_pw92_diff = np.diff(eps_pw92)

# den_min = den_min[(den_min> 0.0001)]
# den_ts = den_ts[(den_ts> 0.0001)]

# t_min = t(trun_rdg_min, trun_den_min)
# t_ts = t(trun_rdg_ts, trun_den_ts)

# A_min = A(den_min)
# A_ts = A(den_ts)

# df_A_min = pd.DataFrame(A_min, columns=['PBE correlation aux. func. GS'])
# df_A_ts = pd.DataFrame(A_ts, columns=['PBE correlation aux. func. TS'])

# fig = plt.figure(dpi=400)

# df_A_min.plot(kind='density',ax = plt.gca())
# df_A_ts.plot(kind='density',ax = plt.gca())

# overlap_index = np.where(abs(A_min - A_ts) < 0.05 )[0]

# A_min[overlap_index] = 0.0
# A_ts[overlap_index] = 0.0

# =============================================================================
# PBE epsilon plot for different t and n
# =============================================================================
# rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

# plt.rcParams.update({'font.size': 12})
    
# fig, ax = plt.subplots(dpi=400)

# cmap = plt.get_cmap('copper',99)
# for i in range(1,200):
    
#     plt.scatter(t_space, t_n_arr_min[:,i], s=2, c=cmap(i), label='Ground state')
#     # plt.scatter(t_space, t_n_arr_ts[:,i], s=2, c=cmap(i), label='Transtion state' )
# norm = mpl.colors.Normalize(vmin=0,vmax=0.05)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# plt.colorbar(sm, ticks=np.linspace(0,0.05,10), 
#                boundaries=np.arange(0.0,0.05,.001))
# plt.plot(t_space, np.sum(0.04*t_n_arr_min, axis=1), color='red', lw=0.5,   label='Ground state' )
# # plt.plot(t_space, np.sum(0.04*t_n_arr_ts, axis=1), color='red', lw=0.5,   label='T state' )
# # ax.set_xlabel(r'Electronic density $\rho(r)$')
# # ax.set_ylabel('$q_0(r)$')
# plt.xlim(0, 3.0)
# plt.ylim(0, 2.5)
# # ax.legend(loc='best', bbox_to_anchor=(0.8, -0.25))
# # ax.legend(loc='best')
# ax.set_xlabel('$s(r)$')
# ax.set_ylabel('Number of grid points')






# =============================================================================
# 
# =============================================================================
s = np.linspace(0, 6, 1000)

# # # s = np.linspace(0, 8, 200)
# # # kappa_PBE=0.804
# # # mu_PBE=0.21951

fig, ax = plt.subplots(dpi=400)

# plt.scatter(trun_den_min, t_min, s=0.05, color='darkviolet', label='G state' )
# plt.scatter(trun_den_ts, t_ts, s=0.05, color='red', label='T state' )
# # plt.ylim(0.0, 4)
# # plt.xlim(0.0, 0.1)
# for n in [0.001, 0.01, 0.05]:
plt.plot(s, gga.eps_c_PBE(0.001, s), '--' , label='n(r) = 0.001' )
plt.plot(s, gga.eps_c_PBE(0.01, s), label='n(r) = 0.01' )
plt.plot(s, gga.eps_c_PBE(0.05, s), label='n(r) = 0.05')
# # plt.plot(n[:-1], eps_pw92_diff, color='darkviolet', label='d\epsilon/dn' )
# # plt.plot(s, gga.F_revPBE_x(s))
# # plt.plot(t, pbe_corr_kernel(t, 0.4))

# # plt.scatter(range(len(A_ts)), A_ts, s=0.1, color='darkviolet', label='T state' )
ax.set_xlabel('$s(r)$')
ax.set_ylabel('$\epsilon^{PBE}_c(n, t)$')
# # ax.set_ylabel('$\epsilon^{PW92} LDA correlation$' )
# # ax.set_ylabel('Gradient part of PBE correlation ($H(t)$)' )
ax.legend(loc='best')