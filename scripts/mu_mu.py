#!/usr/bin/env python3
#  
# COMPUTE_bn.py
# 
# Copyright (c) Rafael L. Delgado
# 

import sys
sys.path.append('..')

import mpi4py
from mpi4py import MPI

import vegas
import numpy as np
from amp import amp

cases={2: (1, 'amp2p', 1000)}

from timeit import default_timer as timer

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ini():
    amp.only_phase_space=False
    amp.indisting_final_state_particles=False
    amp.only_integration=False

    # amp.s=7.e0
    # amp.fp=3.e0

    amp.s=1.e0
    amp.fp=1.e0

def amplitude(n=2):
    start = timer()

    i_dim, ncase, neval = cases[n]
    integrand = getattr(amp, ncase)

    int_amp = vegas.Integrator( i_dim * [[0., 1.]])

    try:
        rank0 = int_amp.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print(f'Running case {ncase} with i_dim={i_dim}')

    try:
        res = int_amp(integrand, neval=neval) #, nint=16, neval=32*1024*1024)
        if rank0:
            print('result = %s    Q = %.2f' % (res, res.Q))
            print(res.summary())
        return res[0].mean
    except:
        if rank0:
            print('FAILURE!!!!')
    finally:
        end = timer()
        if rank0:
            print(f'Case: {ncase}; TIME: {end-start}s')
            print()

# #####################################

def main():
    ini()

    sigma = amplitude(2)

    if rank==0:
        print('Results:')
        print(f'sigma = {sigma}\n')

if __name__ == '__main__':
    main()
