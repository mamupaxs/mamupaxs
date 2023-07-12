tmp = (f1*f2*z12*(f1*f2*z12 - f1 - f2 + 5.d-1))/(2.d0*f1*f2*z12 - f1*z1 - f2*z2) + &
	 (f1*f2*z12*(f1*f2*z12 - f1 - f2 + 5.d-1))/(2.d0*f1*f2*z12 - 2.d0*f1 - 2.d0*f2 + &
	 f1*z1 + f2*z2) + (f1*f3*z13*(f1*f3*z13 - f1 - f3 + 5.d-1))/(2.d0*f1*f3*z13 - &
	 f1*z1 - f3*z3) + (f1*f3*z13*(f1*f3*z13 - f1 - f3 + 5.d-1))/(2.d0*f1*f3*z13 - &
	 2.d0*f1 - 2.d0*f3 + f1*z1 + f3*z3) + ((f1*f2*z12 + f1*f3*z13 - f1)*(f1*f2*z12 + &
	 f1*f3*z13 - f1 - f2 - f3 + 5.d-1))/(-2.d0*f1*f2*z12 - 2.d0*f1*f3*z13 + 2.d0*f1 + &
	 f2*z2 + f3*z3 - 1.d0) + ((f1*f2*z12 + f1*f3*z13 - f1)*(f1*f2*z12 + f1*f3*z13 - &
	 f1 - f2 - f3 + 5.d-1))/(-2.d0*f1*f2*z12 - 2.d0*f1*f3*z13 + 2.d0*f1 + 2.d0*f2 + 2.d0*f3 - f2*z2 - f3*z3 - 1.d0)

select case (idfun)
	case (1)
	M = (2.d0*s) * tmp
	case (2)
	M = (2.d0*s) * tmp*tmp
	case default
	M = 0.d0
end select
