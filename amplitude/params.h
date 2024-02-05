! Introduce here any parameter you want to use inside your integrands.
! These parameters can be modified during execution time via Python instructions.
!
! For instance, if you declare inside this file (params.h):
!
!     double precission :: PARAMETER
!
! you could later use the Python instruction
!
!     amp.PARAMETER = 1.d1
!
! on your Python script to set this parameter to 10 during execution time.
! This PARAMETER variable can be used from inside your integrands (M_2to?.h files).
!
! It is even possible to use the Fortran SELECT CASE statement and an integer
! parameter to select and integrand on execution time.
!
! As an example, we are declaring the widely used floating point variables v and eps, as well as
! an integer parameter for integrating an arbitrary power of function B of formula B.16 on [1].
! Actually, n_param=1 and n_param=2 will be required. See formulas 2.11 and 2.12 of [1].
! 

double precision :: v, eps
integer :: n_param
