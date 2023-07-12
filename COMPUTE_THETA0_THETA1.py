#!/usr/bin/env python3

import mpi4py
from mpi4py import MPI

import vegas
import numpy as np
from amp import amp
cases={2: (1, 'higgs2', 1000), 3: (4, 'higgs3', 4*1024**2), 4: (7, 'higgs4', 1024*1024**2), 5: (10, 'higgs5', 16*1024**2)}

from timeit import default_timer as timer

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ini():
    amp.only_phase_space=False
    amp.indisting_final_state_particles=False
    amp.idfun=1
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
    except:
        if rank0:
            print('FAILURE!!!!')
    finally:
        end = timer()
        if rank0:
            print(f'Case: {ncase}; TIME: {end-start}s')
            print()

        # return res[0].mean
        return res[0].mean

# #####################################

def main():
    ini()

    if rank==0: print('Computing B...')
    amp.idfun=1
    B = amplitude(4)

    if rank==0: print('Computing B^2...')
    amp.idfun=2
    B2 = amplitude(4)
    
    if rank==0: print('Computing phase space...')
    amp.only_phase_space=True
    phas = amplitude(4)

    if rank==0: print(f'B = {B}\nB^2 = {B2}\nphas = {phas}\n')

if __name__ == '__main__':
    main()
