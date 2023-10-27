*************************************************************
## MaMuPaXS
 - Version 1.0, November, 2023
 - Author(s):  Rafael L. Delgado
 - Email:  rafael.delgado@upm.es
*************************************************************

 DESCRIPTION
-------------------------------------------------------------
Evaluating total cross sections from differential cross sections
for 2->n processes, with n>2, can be challenging. This code is
intented for solving such problem for massless initial and
final state particles and n=3,4 and 5. 

Part of the code is Fortran code that compiles to a Python module
via F2PY. For the multidimensional numerical integration, the
VEGAS algorithm is used via the class vegas.Integrator from the
vegas Python module (G.P.Lepage, J.Comput.Phys.**27** (1978) 192
and J.Comput.Phys.**439** (2021) 110386).

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
We have include the Python script *phase_space.py* that computes
the phase space for the 2->n processes with n=2,3,4 and 5. This is
intended to be a test of the code as well as an usage example.

The compilation process generates a Python 3 module on the compiling
directory. Such a module should be loaded into Python via the instruction

```
    from amp import amp
```

Afterwards, the instructions
```
    amp.only_phase_space=True
    amp.v=0.256
    amp.s=1.e-0
```
load v=256MeV, s=1TeV on the evaluation routines. The variable
```only_phase_space``` is set to ```True``` only for evaluating
the phase space. If you are to compute a total cross secction,
```amp.only_phase_space=False``` should be used instead.

If the MPI option is used, the Python script should be called with
the appropriate mpirun command,
```
mpirun -N 6 ./phase_space.py
```
where 6 should be changed to match the number of computing cores
and the Python script should be made executable.

The matrix elements are coded inside the files

* **M_2to2.h**, for a 2->2 process
* **M_2to3.h**, for a 2->3 process
* **M_2to4.h**, for a 2->4 process
* **M_2to5.h**, for a 2->5 process

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

 LICENSE
-------------------------------------------------------------

GNU General Public License (GPLv3).
See detailed text in the [LICENSE](./LICENSE.md) file.

 ATTRIBUTION
-------------------------------------------------------------

We ask that if you use this code for work which results
in a publication that you cite the following paper:

* Rafael L. Delgado, Raquel Gómez-Ambrosio, Javier Martínez-Martín,
Alexandre Salas-Bernárdez, Juan J. Sanz-Cillero,
``SMEFT vs HEFT: multi-Higgs phenomenology'', to appear on arXiv
