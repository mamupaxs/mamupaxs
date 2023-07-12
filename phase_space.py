#!/usr/bin/env python3

PREC=16
#PREC=256

import vegas
import numpy as np
from amp import amp

from timeit import default_timer as timer

def ini():
    amp.only_phase_space=True
    amp.v=0.256
    amp.s=1.e-0


def higgs2():
    start = timer()

    int_higgs2 = vegas.Integrator(1 * [[0., 1.]])
    try:
        rank0 = int_higgs2.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: Higgs2')

    try:
        res_higgs2 = int_higgs2(amp.higgs2, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_higgs2, res_higgs2.Q))
            print(res_higgs2.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: Higgs2; TIME: {end-start}s')
            print()

def higgs3():
    start = timer()

    int_higgs3 = vegas.Integrator(4 * [[0., 1.]])
    try:
        rank0 = int_higgs3.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: Higgs3')

    try:
        res_higgs3 = int_higgs3(amp.higgs3, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_higgs3, res_higgs3.Q))
            print(res_higgs3.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: Higgs3; TIME: {end-start}s')
            print()

def higgs4():
    start = timer()

    int_higgs4 = vegas.Integrator(7 * [[0., 1.]])
    try:
        rank0 = int_higgs4.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: Higgs4')

    try:
        res_higgs4 = int_higgs4(amp.higgs4, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_higgs4, res_higgs4.Q))
            print(res_higgs4.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: Higgs4; TIME: {end-start}s')
            print()

def higgs5():
    start = timer()

    int_higgs5 = vegas.Integrator(10 * [[0., 1.]])
    try:
        rank0 = int_higgs5.mpi_rank == 0
    except:
        rank0 = True

    if rank0:
        print('Running case: Higgs5')

    try:
        res_higgs5 = int_higgs5(amp.higgs5, nitn=16, neval=PREC*1024*1024)
    except:
        if rank0:
            print('FAILURE!!!!')
    else:    
        if rank0:
            print('result = %s    Q = %.2f' % (res_higgs5, res_higgs5.Q))
            print(res_higgs5.summary())
    finally:
        end = timer()
        if rank0:
            print(f'Case: Higgs5; TIME: {end-start}s')
            print()

def main():
    ini()
    for i in range(2):
        for case in [higgs2, higgs3, higgs4, higgs5]:
            case()

if __name__ == '__main__':
    main()
