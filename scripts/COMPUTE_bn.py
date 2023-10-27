#!/usr/bin/env python3

import sys
sys.path.append('..')

import mpi4py
from mpi4py import MPI

import vegas
import numpy as np
from amp import amp
# cases={2: (1, 'amp2p', 1000), 3: (4, 'amp3p', 4*1024**2), 4: (7, 'amp4p', 512*1024**2), 5: (10, 'amp5p', 16*1024**2)}
cases={2: (1, 'amp2p', 1000), 3: (4, 'amp3p', 4*1024), 4: (7, 'amp4p', 512*1024*32), 5: (10, 'amp5p', 16*1024)}

from timeit import default_timer as timer

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ini():
    amp.only_phase_space=False
    amp.indisting_final_state_particles=False
    amp.only_integration=True
    amp.n_param=0
    amp.v=246e-3
    amp.s=1.e0
    amp.eps=0.7
    amp.d=1.23

def amplitude(n=4):
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

    B = {}

    for n in [0, 1, 2, 3, 4]:
        if rank==0: print(f'\n\n----\nComputing B^{n}...')

        amp.n_param = n
        B[n] = amplitude(4)

    if rank==0:
        print('Results:')
        for n,Bn in B.items():
            print(f'B^{n} = {Bn}\n')

if __name__ == '__main__':
    main()
