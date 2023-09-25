#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 13:18:22 2021

@author: mojtaba
"""

from ase.io import read
import os
from ase.calculators.vasp import Vasp
from shutil import copyfile
import numpy as np




Vasp.xc_defaults['pbe'] = {'gga': 'PE', 'lasph': True}
Vasp.xc_defaults['pbe-d3'] = {'gga': 'PE', 'lasph': True, 'ivdw': 11}
Vasp.xc_defaults['pbe-d4'] = {'gga': 'PE', 'lasph': True, 'ivdw': 13}
Vasp.xc_defaults['pbe-ts'] = {'gga': 'PE', 'lasph': True, 'ivdw': 2}
Vasp.xc_defaults['pbe-rev'] = {'gga': 'RE', 'lasph': True}
Vasp.xc_defaults['pbe0'] = {'gga': 'PE', 'lhfcalc': True, 'lasph': True}
Vasp.xc_defaults['vdw-df'] = {'gga': 'RE', 'luse_vdw': True, 'aggac':0.000, 'lasph': True}
Vasp.xc_defaults['vdw-df2'] = {'gga': 'ML', 'luse_vdw': True, 'aggac':0.000, 'lasph': True, 'zab_vdw':-1.8867}
Vasp.xc_defaults['rev-vdw-df2'] = {'gga': 'MK', 'luse_vdw': True, 'aggac':0.000, 'lasph': True, 'zab_vdw':-1.8867, 'param1': 0.1234, 'param2': 0.711357 }
Vasp.xc_defaults['rev-vdw-df2-20'] = {'gga': 'MK', 'luse_vdw': True, 'lhfcalc': True, 'aexx':0.2, 'aggac':0.000, 'lasph': True, 'zab_vdw':-1.8867, 'param1': 0.1234, 'param2': 0.711357}
Vasp.xc_defaults['rev-vdw-df2-25'] = {'gga': 'MK', 'luse_vdw': True, 'lhfcalc': True, 'aexx':0.25, 'aggac':0.000, 'lasph': True, 'zab_vdw':-1.8867, 'param1': 0.1234, 'param2': 0.711357}
Vasp.xc_defaults['cx'] = {'gga': 'CX', 'luse_vdw': True, 'aggac':0.000, 'lasph': True}
Vasp.xc_defaults['cx0'] = {'gga': 'CX', 'luse_vdw': True, 'lhfcalc': True, 'lasph': True, 'aggac':0.000 }
Vasp.xc_defaults['cx0-20'] = {'gga': 'CX', 'luse_vdw': True, 'lhfcalc': True, 'aexx':0.2, 'lasph': True, 'aggac':0.000}
Vasp.xc_defaults['scan'] = {'metagga': 'SCAN', 'lasph': True}
Vasp.xc_defaults['scan-rvv10'] = {'metagga': 'SCAN', 'lasph': True, 'luse_vdw': True, 'bparam': 15.7, 'cparam' : 0.0093}
Vasp.xc_defaults['rvv10'] = {'gga': 'ML', 'lasph': True, 'luse_vdw': True, 'bparam': 6.3, 'cparam' : 0.0093}
Vasp.xc_defaults['b3lyp'] = {'gga': 'B3', 'aexx': 0.2, 'aggac':0.81, 'lhfcalc': True, 'aggax': 0.72, 'aldac': 0.19, 'lasph': True}
Vasp.xc_defaults['optB88-vdW'] = {'gga': 'BO', 'param1': 0.1833333333, 'param2': 0.2200000000, 'luse_vdw': True, 'aggac':0.000, 'lasph': True}
Vasp.xc_defaults['optPBE-vdW'] = {'gga': 'OR', 'luse_vdw': True, 'aggac':0.000, 'lasph': True}
Vasp.xc_defaults['optB86b-vdW'] = {'gga': 'MK', 'param1': 0.1234, 'param2': 1.0000, 'luse_vdw': True, 'aggac':0.000, 'lasph': True}



def get_relax_structure(atoms_file, r_k, functional, out_put_dir, scratch_dir ):
    atoms = read(atoms_file)
    reciprocal = np.sqrt(np.sum(np.array(atoms.cell.reciprocal())**2, axis=1))
    k_pt = r_k*reciprocal + 0.5
    k_pt = np.where(k_pt>1, k_pt.astype(int) , 1)
    calculator = Vasp(kpts=k_pt, algo='Conjugate', nelmin=4,
                      ediff=1E-8, ediffg=-0.01, encut=530, lreal=False,
                      ismear=0, lwave=False, nelm=300, prec='Accurate', sigma=0.05, 
                      ibrion=2, nsw=500, isif=2, xc=functional, ncore=32, lplane=False,
                      directory=out_put_dir)

    atoms.set_calculator(calculator)
    calculator.write_input(atoms)
    copyfile(scratch_dir + '/multi.sh', out_put_dir + '/multi.sh')
    os.chdir(out_put_dir)
    with open('KPOINTS', 'r+') as f:
        with open('a', 'w') as new_file:
            for line in f:
                if 'Monkhorst-Pack' in line:
                    line = line.replace('Monkhorst-Pack', 'Gamma')
                    new_file.write(line)
                else:
                    new_file.write(line)
    os.system('mv a KPOINTS')
    os.system('sbatch multi.sh')



