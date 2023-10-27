 ! Formula B.16 of the article
 ! M = (f1*f2*z12*(f1*f2*z12 - f1 - f2 + 5.d-1))/(2.d0*f1*f2*z12 - f1*z1 - f2*z2) + &
 ! 	 (f1*f2*z12*(f1*f2*z12 - f1 - f2 + 5.d-1))/(2.d0*f1*f2*z12 - 2.d0*f1 - 2.d0*f2 + &
 ! 	 f1*z1 + f2*z2) + (f1*f3*z13*(f1*f3*z13 - f1 - f3 + 5.d-1))/(2.d0*f1*f3*z13 - &
 ! 	 f1*z1 - f3*z3) + (f1*f3*z13*(f1*f3*z13 - f1 - f3 + 5.d-1))/(2.d0*f1*f3*z13 - &
 ! 	 2.d0*f1 - 2.d0*f3 + f1*z1 + f3*z3) + ((f1*f2*z12 + f1*f3*z13 - f1)*(f1*f2*z12 + &
 ! 	 f1*f3*z13 - f1 - f2 - f3 + 5.d-1))/(-2.d0*f1*f2*z12 - 2.d0*f1*f3*z13 + 2.d0*f1 + &
 ! 	 f2*z2 + f3*z3 - 1.d0) + ((f1*f2*z12 + f1*f3*z13 - f1)*(f1*f2*z12 + f1*f3*z13 - &
 ! 	 f1 - f2 - f3 + 5.d-1))/(-2.d0*f1*f2*z12 - 2.d0*f1*f3*z13 + 2.d0*f1 + 2.d0*f2 + 2.d0*f3 - f2*z2 - f3*z3 - 1.d0)
 
 
 ! Formula 2.8 and 2.9 of the article
 t1 = int(dble(f1) * z1)
 t2 = int(dble(f2) * z2)
 t3 = 2 * f1
 t4 = t3 * f2 * z12 - t1 - t2
 t5 = int(dble(f3) * z3)
 t6 = t3 * f3 * z13 - t1 - t5
 t7 = int(dble(f4) * z4)
 t1 = t3 * f4 * z14 - t1 - t7
 t3 = 2 * f2
 t8 = t3 * f3 * z23 - t2 - t5
 t2 = t3 * f4 * z24 - t2 - t7
 t3 = 2 * f3 * f4 * z34 - t5 - t7
 t5 = 1 / t6
 t1 = 1 / t1
 t4 = 1 / t4
 t6 = 1 / t8
 t3 = 1 / t3
 t2 = 1 / t2
 
 M = f1 * f2 * f3 * f4 * ((t4 + t3) * z34 * z12 + (t2 + t5) * z24 * z13 + z23 * z14 * (t1 + t6))
 
 M = (M**n_param)

