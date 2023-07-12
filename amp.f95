module amp
  use ieee_arithmetic

  implicit none
  save

  double precision :: v, d, eps
  double precision :: s
  logical :: only_phase_space, indisting_final_state_particles
  integer :: idfun

  double precision, parameter, private :: PI=4.D0*DATAN(1.D0)
  double precision, private :: tmp

  public  :: higgs2, higgs3, higgs4, higgs5, kin_space_n
  private :: kin_3p, kin_4p, kin_mt4p

contains
  
  subroutine kin_3p(c, s, phi, f, de)
    double precision, intent(in) :: phi(2)
    double precision, intent(inout) :: c(3), s(3)
    double precision, intent(out) :: f(3), de

    double precision pc(2), ps(2), t3

    pc = cos(phi)
    ps = sin(phi)

    t3 = dble(s(1)) * dble(s(2)) * dble(-pc(1) * ps(2) + ps(1) * pc(2)) &
            / (-dble(c(1)) * dble(s(2)) * dble(ps(2)) + dble(c(2)) * dble(s(1)) * dble(ps(1)))

    c(3) = 1.d0/dsqrt(1.d0+t3**2)
    s(3) = dabs(t3)*c(3)

    if (t3.lt.0.d0) then
        c(3) = -c(3)
    end if


    f(1) = -(s(2)) * (ps(2)) * (c(3)) / ((c(1)) * (s(2)) * (ps(2)) - (c(2)) * (s(1)) * (ps(1)) &
            + (c(3)) * (s(1)) * (ps(1)) - (s(2)) * (ps(2)) * (c(3)))

    f(2) = (c(3)) * (s(1)) * (ps(1)) / ((c(1)) * (s(2)) * (ps(2)) - (c(2)) * (s(1)) * (ps(1)) &
            + (c(3)) * (s(1)) * (ps(1)) - (s(2)) * (ps(2)) * (c(3)))

    f(3) = (((c(1)) * (s(2)) * (ps(2)) - (c(2)) * (s(1)) * (ps(1))) / (-(s(1)) * ((c(2)) &
            - (c(3))) * (ps(1)) + (ps(2)) * (s(2)) * ((c(1)) - (c(3)))))

    de = ((s(3) * (-pc(1) * ps(2) + ps(1) * pc(2)) * s(2) + ps(1) * (c(2) * c(3) - 1)) * s(1) - s(2) * ps(2) * (c(1) * c(3) - 1)) * f(3)

    ! de = ((-pc(1) * ps(2) + ps(1) * pc(2)) * (c3 * f(3) + s(3)) * s(2) - ps(1) * (f(3) * c(2) * s3 - c3 * c(2) + 1)) &
    !        * s(1) + s(2) * ps(2) * (f(3) * s3 * c(1) - c3 * c(1) + 1)
  end

  subroutine kin_4p(c, s, phi, f, de)
    double precision, intent(in) :: c(4), s(4), phi(3)
    double precision, intent(out) :: f(4), de

    double precision pc(3), ps(3)
    double precision t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14

    pc = cos(phi)
    ps = sin(phi)

    t1 = pc(2) * ps(3) - pc(3) * ps(2)
    t2 = c(2) * s(4)
    t3 = pc(1) * ps(2) - pc(2) * ps(1)
    t4 = (c(4) - c(1)) * t1
    t5 = t4 * s(3)
    t6 = ps(2) * s(4)
    t7 = t6 * (c(3) - c(1))
    t8 = t3 * (c(4) - c(3))
    t9 = pc(1) * ps(3) - pc(3) * ps(1)
    t10 = s(4) * ps(3)
    t11 = t10 * (c(2) - c(1))
    t12 = (c(4) - c(2)) * t9
    t13 = c(3) - c(2)
    t14 = -(t8 * s(1) + t5 + t7) * s(2) + (t12 * s(1) + t11) * s(3) + s(1) * ps(1) * s(4) * t13
    t8 = (-ps(1) * s(4) * t13 - t12 * s(3) + t8 * s(2)) * s(1)
    t4 = (t4 * s(2) - t11) * s(3) + t8 + t7 * s(2)
    t5 = -(t7 + t5) * s(2) - t8 + t11 * s(3)
    t4 = 0.1D1 / t4
    t7 = 0.1D1 / t5
    t8 = 0.1D1 / t14

    f(1) = (-(c(4) * t1 * s(3) + ps(2) * c(3) * s(4)) * s(2) + t2 * s(3) * ps(3)) * t8
    f(2) = (-(c(4) * t9 * s(3) + ps(1) * c(3) * s(4)) * s(1) + t10 * s(3) * c(1)) * t4
    f(3) = (-(c(4) * t3 * s(2) + t2 * ps(1)) * s(1) + t6 * s(2) * c(1))* t7
    f(4) = ((c(3) * t3 * s(2) - c(2) * s(3) * t9) * s(1) + s(2) * c(1) * s(3) * t1) * t7

    de = -t5

    ! Checks
    ! err(1) = abs( sum(f)-1.d0 )
    ! err(1) = abs( sum(f*c) )
    ! err(1) = abs( f(1)*s(1)*pc(1)+f(2)*s(2)*pc(2)+f(3)*s(3)*pc(3)+f(4)*s(4) )
    ! err(1) = abs( f(1)*s(1)*ps(1)+f(2)*s(2)*ps(2)+f(3)*s(3)*ps(3) )
    ! print *,"max_err = ",maxval(err)
  end

  subroutine kin_mt4p(c, s, phi, SP, f, de)
    double precision, intent(in) :: c(4), s(4), phi(3), SP(4)
    double precision, intent(out) :: f(4), de

    double precision pc(3), ps(3)

    double precision t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
    double precision t11, t12, t13, t14, t15, t16, t17, t18, t19, t20
    double precision t21, t22, t23, t24, t25, t26, t27, t28, t29, t30
    double precision t31, t32, t33, t34, t35, t36, t37

    pc = cos(phi)
    ps = sin(phi)

    t1 = pc(3) * ps(2)
    t2 = -ps(3) * pc(2) + t1
    t3 = SP(1) * c(4) - SP(4)
    t4 = SP(3) * pc(2)
    t5 = -SP(2) * ps(2) + t4
    t6 = c(3) * SP(1) - SP(4)
    t7 = t6 * s(4) - c(3) * SP(2)
    t8 = t5 * c(4)
    t9 = SP(3) * pc(3)
    t10 = -SP(2) * ps(3) + t9
    t11 = c(2) * SP(1) - SP(4)
    t12 = t11 * s(4) - c(2) * SP(2)
    t13 = t10 * c(4)
    t14 = c(3) - c(2)
    t15 = t14 * s(4)
    t16 = ps(1) * pc(2)
    t17 = ps(2) * pc(1) - t16
    t18 = -c(1) + c(3)
    t19 = t18 * s(4)
    t20 = pc(1) * c(3)
    t21 = c(4) * t17
    t22 = (c(4) - c(1)) * t2 * s(3)
    t23 = pc(3) * ps(1)
    t24 = ps(3) * pc(1) - t23
    t25 = c(2) - c(1)
    t26 = t25 * s(4)
    t27 = pc(1) * c(2)
    t28 = c(4) * t24
    t29 = t15 * ps(1)
    t16 = -(-(t20 * s(1) - t19) * ps(2) + (t16 * c(3) + t21) * s(1) -t22) * s(2) + (-(t27 * s(1) - t26) * ps(3) + (t23 * c(2) + t28) *s(1)) * s(3) + t29 * s(1)
    t23 = -SP(2) * ps(1) + SP(3) * pc(1)
    t30 = t23 * c(4)
    t31 = c(1) * SP(1) - SP(4)
    t32 = t31 * s(4) - c(1) * SP(2)
    t33 = (c(4) - c(2)) * t24 * s(3)
    t34 = pc(2) * c(1)
    t35 = c(4) * t2
    t36 = t19 * ps(2)
    t1 = -((-pc(2) * c(3) * s(2) + t15) * ps(1) + (t20 * ps(2) - t21)* s(2) + t33) * s(1) + (-(t34 * s(2) + t26) * ps(3) + (t1 * c(1) - t35) * s(2)) * s(3) + t36 * s(2)
    t21 = (c(4) - c(3)) * t17 * s(2)
    t37 = t26 * ps(3) * s(3)
    t28 = -(-(pc(3) * c(2) * s(3) + t15) * ps(1) + (t27 * ps(3) - t28) * s(3) + t21) * s(1) + (-(pc(3) * c(1) * s(3) + t19) * ps(2) + (t34 * ps(3) + t35) * s(3)) * s(2) + t37
    t22 = t36 - t22
    t21 = (t33 - t21 + t29) * s(1)
    t29 = -t22 * s(2) + t21 + t37

    t28 = 0.1D1 / t28
    t1 = 0.1D1 / t1
    t16 = 0.1D1 / t16
    t29 = 0.1D1 / t29

    f(1) = ((t2 * t3 * s(3) - t4 * c(3) - t7 * ps(2) + t8) * s(2) - (-t12 * ps(3) - t9 * c(2) + t13) * s(3) + t15 * SP(3)) * t16
    f(2) = (-(t24 * t3 * s(3) + t20 * SP(3) + t7 * ps(1) - t30) * s(1) - (-t32 * ps(3) - t9 * c(1) + t13) * s(3) + t19 * SP(3)) * t1
    f(3) = (-(t17 * t3 * s(2) + t12 * ps(1) + t27 * SP(3) - t30) * s(1)- (-t32 * ps(2) - t4 * c(1) + t8) * s(2) + t26 * SP(3)) * t28
    f(4) = -(-(-t24 * t11 * s(3) + t17 * t6 * s(2) - t14 * t23) * s(1) - (-t2 * t31 * s(3) + t18 * t5) * s(2) + s(3) * t25 * t10) * t29

    de = t22 * s(2) - t21 - t37
  end

  subroutine higgs2(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch

    integer :: i
    double precision :: th_c, th_s2, z1, z2, z12, f1, f2, factor, M

    if (nbatch.ne.1) then
       stop "Error: routine 2-higgs requires nbatch=1"
    end if

    do i=1,dim
       th_c = 2.d0*x(1,i) - 1.d0

       factor = 1.d0 / ( 8.d0 * PI ) 

       if (.not.only_phase_space) then
           M = 1.d0

           th_s2 = 1.d0 - th_c**2

           z1 = 1.d0 - th_c
           z2 = 1.d0 + th_c

           f1 = 5.d-1
           f2 = 5.d-1

           z12 = 2.d0

           include 'M_2to2.h' !This file contains the actual amplitude and stores it on M variable

           ans(i) = factor * (5.d-1/s) * M*M
           if (indisting_final_state_particles) then
                   ans(i) = ans(i) / 2.d0
           end if
       else
           ans(i) = factor
       end if


       if ( .not.ieee_is_normal(ans(i)) ) then
          !print *,'.not.ieee_is_normal'
          ans(i)=0.d0
       end if
    end do
  end

  subroutine higgs3(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch
 
    integer :: i
    double precision :: de, f(3), th_c(3), th_s(3), phi(3), f1, f2, f3, f4, z1, z2, z3, z12, z13, z23 
    double precision :: factor, M
  
    if (nbatch.ne.4) then
       stop "Error: routine 3-higgs requires nbatch=4"
    end if

    do i=1,dim

       th_c(1:2) = 2.d0*x(1:2,i)-1.d0
       th_s(1:2) = sqrt(1.d0 - th_c(1:2)**2)
       phi(1:2) = 2.d0*PI*x(3:4,i)
       phi(3) = 0.d0
       
       call kin_3p(th_c, th_s, phi(1:2), f, de)

       if (any(f.lt.0.d0)) then
          ans(i) = 0.d0
       else
          factor = s * product(f) * th_s(3) / ( (2.d0*PI)**2 * 2.d0 * dabs(de) ) 

          if (.not.only_phase_space) then
              M = 1.d0

              f1 = f(1)
              f2 = f(2)
              f3 = f(3)

              z1 = 1.d0 - th_c(1)
              z2 = 1.d0 - th_c(2)
              z3 = 1.d0 - th_c(3)

              z12 = 1.d0 - (th_s(1)*th_s(2)*cos(phi(1)-phi(2)) + th_c(1)*th_c(2))
              z13 = 1.d0 - (th_s(1)*th_s(3)*cos(phi(1)-phi(3)) + th_c(1)*th_c(3))
              z23 = 1.d0 - (th_s(2)*th_s(3)*cos(phi(2)-phi(3)) + th_c(2)*th_c(3))

              include 'M_2to3.h' !This file contains the actual amplitude and stores it on M variable

              ans(i) = factor * (5.d-1/s) * M*M
              if (indisting_final_state_particles) then
                      ans(i) = ans(i) / 6.d0
              end if
          else
              ans(i) = factor
          end if


          if ( .not.ieee_is_normal(ans(i)) ) then
             !print *,'.not.ieee_is_normal'
             ans(i)=0.d0
          end if
       end if

    end do
  end

  subroutine higgs4(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch
 
    integer :: i
    double precision :: de, f(4), th_c(4), th_s(4), phi(3), f1, f2, f3, f4, z1, z2, z3, z4, z12, z13, z14, z23, z24, z34
    double precision :: factor, M
  
    if (nbatch.ne.7) then
       stop "Error: routine 4-higgs requires nbatch=7"
    end if

    do i=1,dim

       th_c(1:4) = 2.d0*x(1:4,i)-1.d0
       th_s = sqrt(1.d0 - th_c**2)
       phi(1:3) = 2.d0*PI*x(5:7,i)
       
       call kin_4p(th_c, th_s, phi, f, de)

       if (any(f.lt.0.d0)) then
          ans(i) = 0.d0
       else

          factor = s**2 * product(f) / ( (2.d0*PI)**4 *dabs(de) ) 

          if (.not.only_phase_space) then
              M = 1.d0

              f1 = f(1)
              f2 = f(2)
              f3 = f(3)
              f4 = f(4)

              z1 = 1.d0 - th_c(1)
              z2 = 1.d0 - th_c(2)
              z3 = 1.d0 - th_c(3)
              z4 = 1.d0 - th_c(4)

              z12 = 1.d0 - (th_s(1)*th_s(2)*cos(phi(1)-phi(2)) + th_c(1)*th_c(2))
              z13 = 1.d0 - (th_s(1)*th_s(3)*cos(phi(1)-phi(3)) + th_c(1)*th_c(3))
              z14 = 1.d0 - (th_s(1)*th_s(4)*cos(phi(1)       ) + th_c(1)*th_c(4))
              z23 = 1.d0 - (th_s(2)*th_s(3)*cos(phi(2)-phi(3)) + th_c(2)*th_c(3))
              z24 = 1.d0 - (th_s(2)*th_s(4)*cos(phi(2)       ) + th_c(2)*th_c(4))
              z34 = 1.d0 - (th_s(3)*th_s(4)*cos(phi(3)       ) + th_c(3)*th_c(4))

              include 'M_2to4.h' !This file contains the actual amplitude and stores it on M variable

              !ans(i) = real(M*conjg(M)*factor)
              ans(i) = factor * (5.d-1/s) * M*M
              if (indisting_final_state_particles) then
                      ans(i) = ans(i) / 24.d0
              end if

              !print *,'s = ',s,'  M = ',M,'  factor = ',factor
              !print *,'f = ',f
              !print *,'zi = ',[z1, z2, z3, z4]
              !print *,'z(12 13 14 23 24 34) = ',[z12, z13, z14, z23, z24, z34]
              !print *,''

          else
              ans(i) = factor
          end if


          if ( .not.ieee_is_normal(ans(i)) ) then
             !print *,'.not.ieee_is_normal'
             ans(i)=0.d0
          end if
       end if

    end do
  end

  subroutine higgs5(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch
 
    integer :: i
    double precision :: de, f(4), f_e(1), th_c(5), th_s(5), phi(3), phi_e(1), SP(4)
    double precision :: f1, f2, f3, f4, f5, z1, z2, z3, z4, z5, z12, z13, z14, z15, z23, z24, z25, z34, z35, z45
    double precision :: factor, M
  
    if (nbatch.ne.10) then
       stop "Error: routine 5-higgs requires nbatch=10"
    end if

    do i=1,dim

       th_c(1:5) = 2.d0*x(1:5,i)-1.d0
       th_s = sqrt(1.d0 - th_c**2)
       phi(1:3) = 2.d0*PI*x(6:8,i)
       phi_e(1) = 2.d0*PI*x(9,i)
       f_e(1) = x(10,i)

       SP(1)=1.d0-f_e(1)
       SP(2)=-f_e(1)*th_s(5)*cos(phi_e(1))
       SP(3)=-f_e(1)*th_s(5)*sin(phi_e(1))
       SP(4)=-f_e(1)*th_c(5)

       call kin_mt4p(th_c(1:4), th_s(1:4), phi(1:3), SP, f, de)

       if (any(f.lt.0.d0).or.(de.eq.0.d0)) then
          ans(i) = 0.d0
       else

          factor = s**3 * product(f) * product(f_e) / ( (2.d0*PI)**6 *dabs(de) ) 

          if (.not.only_phase_space) then
              M = 1.d0

              f1 = f(1)
              f2 = f(2)
              f3 = f(3)
              f4 = f(4)
              f5 = f_e(1)

              z1 = 1.d0 - th_c(1)
              z2 = 1.d0 - th_c(2)
              z3 = 1.d0 - th_c(3)
              z4 = 1.d0 - th_c(4)
              z5 = 1.d0 - th_c(5)

              z12 = 1.d0 - (th_s(1)*th_s(2)*cos(phi(1)-phi(2)  ) + th_c(1)*th_c(2))
              z13 = 1.d0 - (th_s(1)*th_s(3)*cos(phi(1)-phi(3)  ) + th_c(1)*th_c(3))
              z14 = 1.d0 - (th_s(1)*th_s(4)*cos(phi(1)         ) + th_c(1)*th_c(4))
              z15 = 1.d0 - (th_s(1)*th_s(5)*cos(phi(1)-phi_e(1)) + th_c(1)*th_c(5))
              z23 = 1.d0 - (th_s(2)*th_s(3)*cos(phi(2)-phi(3)  ) + th_c(2)*th_c(3))
              z24 = 1.d0 - (th_s(2)*th_s(4)*cos(phi(2)         ) + th_c(2)*th_c(4))
              z25 = 1.d0 - (th_s(2)*th_s(5)*cos(phi(2)-phi_e(1)) + th_c(2)*th_c(5))
              z34 = 1.d0 - (th_s(3)*th_s(4)*cos(phi(3)         ) + th_c(3)*th_c(4))
              z35 = 1.d0 - (th_s(3)*th_s(5)*cos(phi(3)-phi_e(1)) + th_c(3)*th_c(5))
              z45 = 1.d0 - (th_s(4)*th_s(5)*cos(-phi_e(1)      ) + th_c(4)*th_c(5))

              include 'M_2to5.h' !This file contains the actual amplitude and stores it on M variable

              !ans(i) = real(M*conjg(M)*factor)
              ans(i) = factor * (5.d-1/s) * M*M
              if (indisting_final_state_particles) then
                      ans(i) = ans(i) / 120.d0
              end if

              !print *,'s = ',s,'  M = ',M,'  factor = ',factor
              !print *,'f = ',f
              !print *,'zi = ',[z1, z2, z3, z4]
              !print *,'z(12 13 14 23 24 34) = ',[z12, z13, z14, z23, z24, z34]
              !print *,''
          else
              ans(i) = factor
          end if


          if ( .not.ieee_is_normal(ans(i)) ) then
             !print *,'.not.ieee_is_normal'
             ans(i)=0.d0
          end if
       end if

    end do
  end



  subroutine kin_space_n(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch
 
    logical :: t_only_phase_space
    integer :: i, n
    double precision :: de, f(4), phi(3), P(4)
    double precision, allocatable :: th_c(:), th_s(:)
    double precision :: factor, M

    n = 4 + (nbatch - 7)/3

    if ( (n.ne.(3*n-5)) .or. n.lt.4 ) then
       stop "Error: invalid number of arguments"
    end if

    allocate(th_c(n),stat=i)
    if (i.ne.0) then
       stop "Error: couldn't allocate memory for array!"
    end if

    allocate(th_s(n),stat=i)
    if (i.ne.0) then
       stop "Error: couldn't allocate memory for array!"
    end if
  
    only_phase_space = only_phase_space
    only_phase_space = .true.

    do i=1,dim

       th_c(1:4) = 2.d0*x(1:4,i)-1.d0
       phi(1:3) = 2.d0*PI*x(n+1:n+3,i)
       if (n.gt.4) then
           ! th_c(5:n) = 2.d0*x(8:3+n,i)-1.d0
           ! phi(4:n) = 2.d0*PI*x(n+1:n+3,i) 
       end if

       th_s = sqrt(1.d0 - th_c**2)

       if (n.eq.4) then
           call kin_4p(th_c, th_s, phi, f, de)
       else if (n.gt.4) then
           ! P(1) = sum()
           call kin_mt4p(th_c(1:4), th_s(1:4), phi(1:3), P, f, de)
       else
           only_phase_space = t_only_phase_space
           deallocate(th_c)
           deallocate(th_s)
           stop "Error: wrong number of arguments!!"
       end if


       if (any(f.lt.0.d0)) then
          ans(i) = 0.d0
       else
          factor = s**2 * product(f) / ( (2.d0*PI)**4 *dabs(de) ) 
          ans(i) = factor

          if ( .not.ieee_is_normal(ans(i)) ) then
             !print *,'.not.ieee_is_normal'
             ans(i)=0.d0
          end if
       end if

    end do
    only_phase_space = t_only_phase_space
    deallocate(th_c)
    deallocate(th_s)
  end
end module amp
