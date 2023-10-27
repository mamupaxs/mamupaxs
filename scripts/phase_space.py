#!/usr/bin/env python3

import sys
sys.path.append('..')

PREC=16
#PREC=256

import vegas
import numpy as np
from amp import amp

from timeit import default_timer as timer

def ini():
    amp.only_phase_space=True
    amp.only_integration=False
    amp.indistinguisable_particles=False
    amp.v=0.256
    amp.s=1.e-0


def amp2p():
    start = timer()

    int_amp2p = vegas.Integrator(1 * [[0., 1.]])
    try:
        rank0 = int_amp2p.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: 2 particles')

    try:
        res_amp2p = int_amp2p(amp.amp2p, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_amp2p, res_amp2p.Q))
            print(res_amp2p.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: 2 particles; TIME: {end-start}s')
            print()

def amp3p():
    start = timer()

    int_amp3p = vegas.Integrator(4 * [[0., 1.]])
    try:
        rank0 = int_amp3p.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: 3 particles')

    try:
        res_amp3p = int_amp3p(amp.amp3p, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_amp3p, res_amp3p.Q))
            print(res_amp3p.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: 3 particles; TIME: {end-start}s')
            print()

def amp4p():
    start = timer()

    int_amp4p = vegas.Integrator(7 * [[0., 1.]])
    try:
        rank0 = int_amp4p.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: 4 particles')

    try:
        res_amp4p = int_amp4p(amp.amp4p, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_amp4p, res_amp4p.Q))
            print(res_amp4p.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: 4 particles; TIME: {end-start}s')
            print()

def amp5p():
    start = timer()

    int_amp5p = vegas.Integrator(10 * [[0., 1.]])
    try:
        rank0 = int_amp5p.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: 5 particles')

    try:
        res_amp5p = int_amp5p(amp.amp5p, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_amp5p, res_amp5p.Q))
            print(res_amp5p.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: 5 particles; TIME: {end-start}s')
            print()

def main():
    ini()
    for i in range(2):
        for case in [amp2p, amp3p, amp4p, amp5p]:
            case()

if __name__ == '__main__':
    main()
