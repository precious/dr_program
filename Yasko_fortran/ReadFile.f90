include 'Lapl_3D.f90'
!------------------------------------------------------------------------!
!-----  (c) Mykola Yas'ko. All right reserved.                    -------!
!-----  http://nyasko.homepage.com   E-mail: nyasko@hotmail.com   -------!
!-----  This is Demo-program for solving 3D Laplace equation      -------!
!-----  by the direct boundary element method.                    -------!
!-----  Domain is interior of cube                                -------!
!------------------------------------------------------------------------!
PROGRAM TEST_BEM       ! Interior problem
  USE BEM_Lapl_3D      ! Analytical solution  \phi = x**2 - y**2
  IMPLICIT NONE
  integer :: N         ! Number of boundary elements. 
  type(Vector3D)     :: M, Grad
  integer :: i, Method=2
  real(bem) :: Pot
  real :: time0, time1

  call CPU_TIME(time0)     ! Remember start time

  N = InputBND3('Cube2352.bnd', .false.)  ! Input boundary elements vertices 1 2 3
                                          ! if .true.   then 1 3 2
  if(N<=0) then
     print*, ' Error=', N, ' Download file Cube2352.bnd'
     stop
  end if

  i = FillBE()     ! Additional computation for boundary elements

  if ( i == -7 ) then
     stop ' Bad boundary elements occured. Stop!'
  elseif ( i == -6 ) then
     stop ' Function InitialiseBE did not used. Stop!'
  end if

!-------- Boundary condition of 1st kind  --- ----------------------
! BE%bc = 1   ! Kind of boundary condition
! BE%POT = BE%P%x**2-BE%P%y**2   ! Analytical solution
!-------------------------------------------------------------------

!-------- Boundary condition of 2nd kind  --- ----------------------
! BE%bc = 2   ! Kind of boundary condition
! BE%Vn = 2*(BE%P%x*BE%n%x-BE%P%y*BE%n%y)     ! Normal derivative
!-------------------------------------------------------------------

!-------- Boundary condition of 1st kind: Potential + Tangent velocity ---
  BE%bc = 5   ! Kind of boundary condition
  BE%POT = BE%P%x**2-BE%P%y**2   ! Analytical solution
  BE%Vt =  BE%n.vp.( (/ (Vector3D(2*BE(i)%P%x,-2*BE(i)%P%y,0._bem), i=1,N) /) .vp. BE%n)
!-------------------------------------------------------------------------

  print*, 'Interior boundary-valued problem with', size(BE),' boundary elements'
  i = SolveBVP(Method)         ! Solve boundary-valued problem by 
  if (i < 0 ) then             ! conjugate gradient method (Method=1)
    select case (i)
      case (-1)
        print*, ' Linear system was not solved'
      case (-10)
        print*, ' Unknown boundary condition'
      case (-5)
        print*, ' Function FillBE did not used'
      case (-4)
        print*, ' Boundary Elements are out of memory'
      case default
        print*, 'Unknown error ', i
    end select
    stop
  end if
  print*, '   Iterations=', i

  print*, '   Max error =', MaxVal(abs(BE%POT-(BE%P%x**2-BE%P%y**2))), &
                            MaxVal(abs(BE%Vn-(2*(BE%P%x*BE%n%x-BE%P%y*BE%n%y) )))

  M = Vector3D( 0.1_bem, 0.2_bem, 0.3_bem)
  call Result( M, Pot, Grad)
  Pot = Potential(M)
  print *, ' Point ', M
  print *, ' Numerical Potential ', Pot
  print *, ' AnalyticalPotential ', M%x**2-M%y**2
  print *, ' Numerical  Gradient ', Gradient(M)
  print *, ' Analytical Gradient ',Vector3D(2*M%x, -2*M%y, 0.0_bem)

  M = Vector3D( 0.2_bem, 0.3_bem, 0.4_bem)
  call Result( M, Pot, Grad)
  Pot = Potential(M)
  print *, ' Point ', M
  print *, ' Numerical Potential ', Pot
  print *, ' AnalyticalPotential ', M%x**2-M%y**2
  print *, ' Numerical  Gradient ', Gradient(M)
  print *, ' Analytical Gradient ',Vector3D(2*M%x, -2*M%y, 0.0_bem)

  call CPU_TIME(time1)
  print*, ' Total computational time was ', time1-time0,'s.'
  print*, '        Normal derivative on boundary'
  print*, '         N         Numer.           Anal.'
   do i=1,N
      print*, i, BE(i)%Vn, 2*(BE(i)%P%x*BE(i)%n%x-BE(i)%P%y*BE(i)%n%y)
   end do
  stop
END PROGRAM TEST_BEM
