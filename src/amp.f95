!  
! amp.f95
! 
! Copyright (c) Rafael L. Delgado
!
! Distributed under the GPL v3 license
!
! Sections of the code (marked accordingly) has been generated
! with SageMath and Sympy, both open source programs.
! 

module amp
  use ieee_arithmetic

  implicit none
  save

  double precision, public :: s
  logical, public :: only_phase_space, only_integration, indisting_final_state_particles


  public  :: amp2p, amp3p, amp4p, amp5p, kin_space_np
  include "../amplitude/params.h"

  private :: kin_3p, kin_4p, kin_mt4p
  double precision, parameter, private :: PI=4.D0*DATAN(1.D0)
  include "../amplitude/tmp_vars.h"

contains
  
  subroutine kin_3p(c, s, phi, f, de)
    double precision, intent(in) :: phi(2)
    double precision, intent(inout) :: c(3), s(3)
    double precision, intent(out) :: f(3), de

    double precision pc(2), ps(2), t3

    double precision x0, x1, x2, x3, x4, x5
    double precision dsA, t1, t2

    pc = cos(phi)
    ps = sin(phi)

    ! section generated with SageMath and Sympy
    x0 = ps(2)*s(2)
    x1 = x0*c(1) - c(2)*ps(1)*s(1)
    x2 = -x1
    x3 = (pc(1)*ps(2) - pc(2)*ps(1))*s(1)*s(2)
    t1 = x2
    t2 = x3
    ! end of section

    s(3) = t2/sqrt(t1**2+t2**2)
    c(3) = -t1/sqrt(t1**2+t2**2)

    if (s(3).lt.0.d0) then
        s(3) = -s(3)
        c(3) = -c(3)
    end if

    ! section generated with SageMath and Sympy
    x4 = ps(1)*s(1)
    x5 = -x0 + x4
    dsA = -x1 - x5*c(3)
    f(1) = x0*c(3)
    f(2) = -x4*c(3)
    f(3) = x2
    ! end of section
    
    f = f / dsa
    
    ! section generated with SageMath and Sympy
    de = f(3)*(-x1*c(3) - x3*s(3) - x5)
    ! end of section
  end

  subroutine kin_4p(c, s, phi, f, de)
    double precision, intent(in) :: c(4), s(4), phi(3)
    double precision, intent(out) :: f(4), de

    double precision pc(3), ps(3)
    double precision x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12

    pc = cos(phi)
    ps = sin(phi)

    ! section generated with SageMath and Sympy
    x0 = pc(1)*ps(2) - pc(2)*ps(1)
    x1 = s(1)*s(2)
    x2 = x0*x1
    x3 = x2*c(3)
    x4 = pc(1)*ps(3) - pc(3)*ps(1)
    x5 = x4*s(1)
    x6 = pc(2)*ps(3) - pc(3)*ps(2)
    x7 = (-x5*c(2) + x6*c(1)*s(2))*s(3)
    x8 = ps(1)*s(1)
    x9 = ps(2)*s(2)
    x10 = -x8*c(2) + x9*c(1)
    x11 = c(4)*s(3)
    x12 = ps(3)*s(3)
    de = -x3 - x7 + (x0*x1 + (-x4*s(1) + x6*s(2))*s(3))*c(4) + (-x10 - (x8 &
          - x9)*c(3) + (c(1) - c(2))*ps(3)*s(3))*s(4)
    f(1) = x11*x6*s(2) - (x12*c(2) - x9*c(3))*s(4)
    f(2) = -x11*x5 + (x12*c(1) - x8*c(3))*s(4)
    f(3) = -x10*s(4) + x2*c(4)
    f(4) = -x3 - x7
    ! end of section

    f = f / de
  end

  subroutine kin_mt4p(c, s, phi, SP, f, de)
    double precision, intent(in) :: c(4), s(4), phi(3), SP(4)
    double precision, intent(out) :: f(4), de

    double precision pc(3), ps(3)

    double precision x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
    double precision x11, x12, x13, x14, x15, x16, x17, x18, x19, x20
    double precision x21, x22, x23, x24, x25, x26, x27

    pc = cos(phi)
    ps = sin(phi)

    ! section non-usage for commercial purposes
    x0 = pc(1)*ps(2) - pc(2)*ps(1)
    x1 = pc(1)*ps(3) - pc(3)*ps(1)
    x2 = x1*s(1)
    x3 = pc(2)*ps(3)
    x4 = x3 - pc(3)*ps(2)
    x5 = ps(2)*s(2)
    x6 = ps(1)*s(1)
    x7 = -x0*s(1)
    x8 = -x1*s(1)
    x9 = -x4
    x10 = x9*s(2)
    x11 = SP(2)*ps(2) - SP(3)*pc(2)
    x12 = x11*s(2)
    x13 = SP(2)*ps(3) - SP(3)*pc(3)
    x14 = -SP(4)
    x15 = ps(3)*s(3)
    x16 = -SP(3)
    x17 = SP(2)*ps(1) - SP(3)*pc(1)
    x18 = x13*c(1) - x2*SP(4)
    x19 = x16 + x6*SP(1)
    x20 = SP(1)*c(1)
    x21 = x14 + x20
    x22 = x6*SP(4) - SP(3)*c(1)
    x23 = x17*s(1)
    x24 = x13 + x8*SP(1)
    x25 = -x0*SP(4)*s(1) + x11*c(1)
    x26 = x23 - (x11 + x7*SP(1))*s(2)
    de  = -x0*c(3)*s(1)*s(2) - (-x2*c(2) + x4*c(1)*s(2))*s(3) + (-x7*s(2) + &
          (-x10 + x8)*s(3))*c(4) + (-x5*c(1) - (-x5 + x6)*c(3) + (c(1) - c( &
          2))*ps(3)*s(3) + c(2)*ps(1)*s(1))*s(4)
    f(1) = -x12*c(3) - (-x12 + (x10*SP(1) + x13)*s(3))*c(4) - (-x13*c(2) + &
          x4*SP(4)*s(2))*s(3) - (x15*(x14 + SP(1)*c(2)) + x5*SP(4) - (x16 + &
          x5*SP(1))*c(3) - SP(3)*c(2))*s(4)
    f(2) = x17*c(3)*s(1) - x18*s(3) - (x23 - x24*s(3))*c(4) + (x15*x21 - x19 &
          *c(3) + x22)*s(4)
    f(3) = -x23*c(2) + x25*s(2) + x26*c(4) - (-x19*c(2) + x21*x5 + x22)*s(4)
    f(4) = x17*c(2)*s(1) - x25*s(2) - x26*c(3) + (x18 - x24*c(2) + (x20*x9 + &
          x3*SP(4) - SP(4)*pc(3)*ps(2))*s(2))*s(3)
    ! end of section non-usable for commercial purposes

    f = f / de
  end

  subroutine amp2p(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch

    integer :: i
    double precision :: th_c, th_s2, z1, z2, z12, f1, f2, factor, M

    if (nbatch.ne.1) then
       stop "Error: routine 2-particle requires nbatch=1"
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

           include "../amplitude/M_2to2.h" !This file contains the actual amplitude and stores it on M variable

           if (only_integration) then
                ans(i) = M*factor
           else
                ans(i) = factor * (5.d-1/s) * M*M
                if (indisting_final_state_particles) then
                        ans(i) = ans(i) / 2.d0
                end if
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

  subroutine amp3p(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch
 
    integer :: i
    double precision :: de, f(3), th_c(3), th_s(3), phi(2), f1, f2, f3, f4, z1, z2, z3, z12, z13, z23 
    double precision :: factor, M
  
    if (nbatch.ne.4) then
       stop "Error: routine 3-particle requires nbatch=4"
    end if

    do i=1,dim

       th_c(1:2) = 2.d0*x(1:2,i)-1.d0
       th_s(1:2) = sqrt(1.d0 - th_c(1:2)**2)
       phi(1:2) = 2.d0*PI*x(3:4,i)
       
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
              z13 = 1.d0 - (th_s(1)*th_s(3)*cos(phi(1)       ) + th_c(1)*th_c(3))
              z23 = 1.d0 - (th_s(2)*th_s(3)*cos(phi(2)       ) + th_c(2)*th_c(3))

              include "../amplitude/M_2to3.h" !This file contains the actual amplitude and stores it on M variable

              if (only_integration) then
                   ans(i) = real(M)*factor
              else
                   ans(i) = factor * (5.d-1/s) * M*M
                   if (indisting_final_state_particles) then
                           ans(i) = ans(i) / 6.d0
                   end if
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

  subroutine amp4p(ans, x, dim, nbatch)
    integer, intent(in) :: dim, nbatch
    double precision, intent(in) :: x(nbatch, dim)
    double precision, intent(out) :: ans(dim)
    ! f2py intent(out) :: ans
    ! f2py intent(in) :: x, dim, nbatch
 
    integer :: i
    double precision :: de, f(4), th_c(4), th_s(4), phi(3), f1, f2, f3, f4, z1, z2, z3, z4, z12, z13, z14, z23, z24, z34
    double precision :: factor, M
  
    if (nbatch.ne.7) then
       stop "Error: routine 4-particle requires nbatch=7"
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

              include "../amplitude/M_2to4.h" !This file contains the actual amplitude and stores it on M variable

              if (only_integration) then
                   ans(i) = real(M)*factor
              else
                   ans(i) = factor * (5.d-1/s) * M*M
                   if (indisting_final_state_particles) then
                           ans(i) = ans(i) / 24.d0
                   end if
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

  subroutine amp5p(ans, x, dim, nbatch)
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
       stop "Error: routine 5-particle requires nbatch=10"
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

              include "../amplitude/M_2to5.h" !This file contains the actual amplitude and stores it on M variable

              if (only_integration) then
                   ans(i) = real(M)*factor
              else
                   ans(i) = factor * (5.d-1/s) * M*M
                   if (indisting_final_state_particles) then
                           ans(i) = ans(i) / 120.d0
                   end if
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



  subroutine kin_space_np(ans, x, dim, nbatch)
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
