! Introduce your 2->4 matrix element or integrand here.
! For illustrative purposes, we have set the identity function
! as integrand.
! 
! M = 1.d0
!
! The previous definition of the integrand M has been commented out
! in order to introduce formula B.16 of article [1]:

x0 = f1*z1
x1 = f2*z2
x2 = 2*f1
x3 = z12*z34
x4 = f3*z3
x5 = f4*z4
x6 = 2*f4
x7 = z13*z24
x8 = z14*z23
x9 = f2*f3
M = -f1*f4*x9*(x3/(-f3*x6*z34 + x4 + x5) + x3/(-f2*x2*z12 + x0 + x1) &
      + x7/(-f3*x2*z13 + x0 + x4) + x7/(-f2*x6*z24 + x1 + x5) + x8/(-f4 &
      *x2*z14 + x0 + x5) + x8/(x1 + x4 - 2*x9*z23))

M = M**n_param

