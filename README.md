*************************************************************
## MaMuPaXS
 - Version 1.1, November, 2023
 - Author(s):  Rafael L. Delgado
 - Email:  rafael.delgado@upm.es
*************************************************************

 DESCRIPTION
-------------------------------------------------------------
Evaluating total cross sections from differential cross sections
for 2->n processes, with n>2, can be challenging. This code is
intended for solving such problem for massless initial and
final state particles and n=3,4 and 5. Furthermore, the framework
is kept as simple as possible for this single task.

Part of the code is Fortran code that compiles to a Python module
via F2PY. For the multidimensional numerical integration, the
VEGAS algorithm is used via the class vegas.Integrator from the
vegas Python module (G.P.Lepage, J.Comput.Phys.**27** (1978) 192
and J.Comput.Phys.**439** (2021) 110386).

The code is described on Ref. [1]. Furthermore, another
upcoming article will describe it in detail.

 REQUIREMENTS
-------------------------------------------------------------
This code requires a Fortran compiler and a Python 3 installation.
The Python packages vegas, Numpy and F2PY should be available.
Optionally, MPI can be used for parallelizing the Monte Carlo
integration. For using this option, the mpi4py Python module
is also required.

 COMPILING
-------------------------------------------------------------
The repository includes a Makefile. The only required command
is:

```make```

 USAGE
-------------------------------------------------------------
The repository contains the following folders:
* **scripts**: two Python scripts as an example of usage
* **amplitude**: Fortran source files with the amplitudes coded by the user
* **src**: the actual source code

The compilation process generates a Python 3 module on the repository
directory. Such a module should be loaded into Python via the instruction

```
    from amp import amp
```

Afterwards, the instructions
```
    amp.only_phase_space=False
    amp.indisting_final_state_particles=False
    amp.only_integration=False

    amp.s=1.e-0
```
load s=1TeV on the evaluation routines. There are three boolean varibles
for special purposes:
* **amp.only_phase_space**: if set to True, only the phase space is integrated
* **amp.indisting_final_state_particles**: if set to True, the amplitudes
* will be divided by $n!$. This is required if all the final state
* particles are indistinguishable.
* **amp.only_integration**: instructs MaMuPaX to directly integrate the value
* of the amplitude instead of treating it as a quantum field theory amplitude.
* The latter requires to compute the square modulus and divide by $2s$.

If the MPI option is used, the Python script should be called with
the appropriate mpirun command,
```
mpirun -N 6 ./phase_space.py
```
where 6 should be changed to match the number of computing cores. Note that
the Python script should be made executable.

The folder **scripts**, contains the Python script *phase_space.py*
that computes the phase space for the 2->n processes with n=2,3,4 and 5. This is
intended to be a test of the code as well as an usage example. It also contains
the script *COMPUTE_bn.py*, for computing the integration of the $B^n$ variable
($n=0,1,2,3,4$) of equations (2.8) and (2.9) in our reference [1].

The matrix elements should be coded in Fortran by the user inside the following files:

* **amplitude/M_2to2.h**, for a 2->2 process
* **amplitude/M_2to3.h**, for a 2->3 process
* **amplitude/M_2to4.h**, for a 2->4 process
* **amplitude/M_2to5.h**, for a 2->5 process

In each of these files, the code should return a value with
the matrix element on the ```M``` variable. The following variables
describing the outgoing particles 4-momenta can be used (notation: i,j=1,2,3,4,5):

* **fi**      : $\lVert\vec p_i\rVert/\sqrt{s}$, three-momentum fractions
* **zi**      : $1-\cos\theta_i$, angular functions
* **zij**     : $1-\cos\theta_{ij}$, where $\theta_{ij}$ is the angle between
  the $i$-th and $j$-th outgoing particle

The following variables can also be used, although were originally intended for
MaMuPaXS internal use:
* **th_c(i)** : $\cos\theta_i$, where $\theta_i$ is the azimuthal angle of
  the $i$-th outgoing particle 
* **th_s(i)** : $\sin\theta_i$
* **phi(i)**  : $\phi_i$, polar angle ($i\leq 3$)
* **phi_e(k)**: $\phi_i$, polar angle for $k=i-4$, where $i$ refers to
  the $i$-th outgoing particle ($i\geq 5$).

Note that the code assumes, without loss of generality, that $\phi_l=0$, where
$l=3$ for 3 outgoing particles and $l=4$ for $\geq$ 4 outgoing particles.

Furthermore, the user can declare parameters and temporary variables in these two files:
* **amplitude/params.h**: parameters that can be configured inside the Python code and
  accesed by the user's Fortran code that encodes the amplitude.
* **amplitude/tmp_vars.h**: temporary variables used by the user's Fotran code.

 LICENSE
-------------------------------------------------------------

In this version of the program, the sections that were generated by Maple software
have been re-generated with SageMath [2] and SymPy [3]. Since SageMath and SymPy
are available under the GPL and BSD licenses, this allows for the distribution of
this version of MaMuPaXS under a GPL license (see detailed text in the [LICENSE.md](LICENSE.md)
file).

Please, note that the old version 1.0 of MaMuPaXS was not available under the GPL license
due to part of the analytical expressions being generated by Maple software,
hence making it unusable for commercial purposes according to:
https://www.maplesoft.com/documentation_center/Maplesoft_EULA.pdf

 ATTRIBUTION
-------------------------------------------------------------

We ask that if you use this code for work which results
in a publication that you cite the following paper:

* [1] Rafael L. Delgado, Raquel Gómez-Ambrosio, Javier Martínez-Martín,
Alexandre Salas-Bernárdez, Juan J. Sanz-Cillero,
``SMEFT vs HEFT: multi-Higgs phenomenology'', to appear on arXiv


 BIBLIOGRAPHY
-------------------------------------------------------------
For the generation of several sections of the code that compute analytical expressions
we have used SageMath and the fcode function of Sympy. See our publication [1] for details
and the following referentes for SageMath and Sympy:


* [2] SageMath, the Sage Mathematics Software System (Version 10.1),
   The Sage Developers, 2023, https://www.sagemath.org. 

* [3] Meurer A, Smith CP, Paprocki M, Čertík O, Kirpichev SB, Rocklin M, Kumar A,
Ivanov S, Moore JK, Singh S, Rathnayake T, Vig S, Granger BE, Muller RP,
Bonazzi F, Gupta H, Vats S, Johansson F, Pedregosa F, Curry MJ, Terrel AR,
Roučka Š, Saboo A, Fernando I, Kulal S, Cimrman R, Scopatz A. (2017) SymPy:
symbolic computing in Python. *PeerJ Computer Science* 3:e103
https://doi.org/10.7717/peerj-cs.103
