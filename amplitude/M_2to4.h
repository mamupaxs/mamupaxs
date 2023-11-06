 ! Formula B.16 of the article
  M = (f1*f2*z12*(f1*f2*z12 - f1 - f2 + 5.d-1))/(2.d0*f1*f2*z12 - f1*z1 - f2*z2) + &
  	 (f1*f2*z12*(f1*f2*z12 - f1 - f2 + 5.d-1))/(2.d0*f1*f2*z12 - 2.d0*f1 - 2.d0*f2 + &
  	 f1*z1 + f2*z2) + (f1*f3*z13*(f1*f3*z13 - f1 - f3 + 5.d-1))/(2.d0*f1*f3*z13 - &
  	 f1*z1 - f3*z3) + (f1*f3*z13*(f1*f3*z13 - f1 - f3 + 5.d-1))/(2.d0*f1*f3*z13 - &
  	 2.d0*f1 - 2.d0*f3 + f1*z1 + f3*z3) + ((f1*f2*z12 + f1*f3*z13 - f1)*(f1*f2*z12 + &
  	 f1*f3*z13 - f1 - f2 - f3 + 5.d-1))/(-2.d0*f1*f2*z12 - 2.d0*f1*f3*z13 + 2.d0*f1 + &
  	 f2*z2 + f3*z3 - 1.d0) + ((f1*f2*z12 + f1*f3*z13 - f1)*(f1*f2*z12 + f1*f3*z13 - &
  	 f1 - f2 - f3 + 5.d-1))/(-2.d0*f1*f2*z12 - 2.d0*f1*f3*z13 + 2.d0*f1 + 2.d0*f2 + 2.d0*f3 - f2*z2 - f3*z3 - 1.d0)
 
 
 ! Formula 2.8 and 2.9 of the article
 !M = f1 * f2 * f3 * f4 * (z12 * z34 / (2 * f1 * f2 * z12 - f1 * z1 - f2 * z2) + z13 * z24 / (2 * f1 * f3 * z13 - f1 * z1 - f3 * z3) + z14 * z23 / (2 * f1 * f4 * z14 - f1 * z1 - f4 * z4) + z23 * z14 / (2 * f2 * f3 * z23 - f2 * z2 - f3 * z3) + z24 * z13 / (2 * f2 * f4 * z24 - f2 * z2 - f4 * z4) + z34 * z12 / (2 * f3 * f4 * z34 - f3 * z3 - f4 * z4))

 
 M = M**n_param

