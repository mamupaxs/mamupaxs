*************************************************************
## MaMuPaXS
 - Version 1.0, July, 2023
 - Author(s):  Rafael L. Delgado
 - Email:  rafael.delgado@upm.es
*************************************************************

 DESCRIPTION
-------------------------------------------------------------
Evaluating total cross sections from differential cross sections
for 2->n processes, with n>2, can be challenging. This code is
intented for solving such problem for massless initial and
final state particles and n=3,4 and 5. In 

Part of the code is Fortran code that compiles to a Python module
via F2PY. For the multidimensional numerical integration, the
VEGAS algorithm is used via the class vegas.Integrator from the
vegas Python module (G.P.Lepage, J.Comput.Phys.**27** (1978) 192
and J.Comput.Phys.**439** (2021) 110386).

 REQUIREMENTS
-------------------------------------------------------------
This code requires a working Fortran compiler and a Python 3
installation. The Python packages vegas, Numpy and F2PY should
be available.

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

 LICENSE
-------------------------------------------------------------

GNU General Public License (GPLv3).
See detailed text in the [LICENSE](./LICENSE) file.

 ATTRIBUTION
-------------------------------------------------------------
