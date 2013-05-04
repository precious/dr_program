include 'Lapl_3D.f90'
!------------------------------------------------------------------------!
!-----  (c) Mykola Yas'ko. All right reserved.                    -------!
!-----  http://nyasko.homepage.com   E-mail: nyasko@hotmail.com   -------!
!-----  This is Demo-program for solving 3D Laplace equation      -------!
!-----  by the direct boundary element method.                    -------!
!-----  Domain is interior of sphere, radius=1                    -------!
!------------------------------------------------------------------------!
PROGRAM TEST_BEM       ! HemiSphere, 1024 BE, interior problem
  USE BEM_Lapl_3D      ! Analytical solution  \phi = z
  IMPLICIT NONE        ! Symmetry plane z=0 Pot=0
  integer, parameter :: N = 1024 ! Number of boundary elements. 
  type(Vector3D)     :: M, Grad
  integer :: i, Method=4
  real(bem) :: Pot
  real :: time0, time1

  call CPU_TIME(time0)     ! Remember start time

!---------- Allocate N boundary elements ---------
  i=InitialiseBE(N)
!--------------------------------------------------------

!-------- Input boundary elements vertices --------------
  open(7,file="s2048.dat", position="rewind", status='old', iostat=i)
  if(i/=0) stop ' Bad or missing data file: s2048.dat'
  do i=1,N
     read(7,*) BE(i)%P1, BE(i)%P2, BE(i)%P3 ! Interior problem in this order
  end do
  close(7)
!--------------------------------------------------------

  i = FillBE()     ! Additional computation for boundary elements

  if ( i == -7 ) then
     stop ' Bad boundary elements occured. Stop!'
  elseif ( i == -6 ) then
     stop ' Function InitialiseBE did not used. Stop!'
  end if

  i = SetSymmetry( 0, 0._bem,  0, 0._bem, -1, 0._bem) ! Plane symmetry z=0, Pot=0

!-------- Boundary condition of 2nd kind  --- ----------------------
  BE%bc = 2       ! Kind of boundary condition
  BE%Vn = (/ (Vector3D(BE(i)%P%y*BE(i)%P%z, BE(i)%P%x*BE(i)%P%z, BE(i)%P%x*BE(i)%P%y),i=1,N) /).sp.BE%n  ! Normal derivative
!-------------------------------------------------------------------

!-------- Boundary condition of 1nd kind  --- ----------------------
! BE%bc = 1       ! Kind of boundary condition
! BE%POT = BE%P%z ! Analytical solution
!-------------------------------------------------------------------

  print*, 'Interior of hemisphere with',size(BE),' boundary elements'
  i = SolveBVP(Method)         ! Solve boundary-valued problem by 
  if (i < 0 ) then             ! conjugate gradient method
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

  print*, '   Max error =', MaxVal(abs(BE%POT-BE%P%z)), &
                            MaxVal(abs(BE%Vn-BE%n%z))

  M = Vector3D( 0.1_bem, 0.2_bem, 0.3_bem)
  call Result( M, Pot, Grad)
  Pot = Potential(M)
  print *, ' Point ', M
  print *, ' Numerical Potential ', Pot
  print *, ' AnalyticalPotential ', M%x*M%y*M%z
  print *, ' Numerical  Gradient ', Gradient(M)
  print *, ' Analytical Gradient ',Vector3D( M%y*M%z, M%x*M%z, M%x*M%y)

  M = Vector3D( 0.2_bem, 0.3_bem, 0.1_bem)
  call Result( M, Pot, Grad)
  Pot = Potential(M)
  print *, ' Point ', M
  print *, ' Numerical Potential ', Pot
  print *, ' AnalyticalPotential ', M%x*M%y*M%z
  print *, ' Numerical  Gradient ', Gradient(M)
  print *, ' Analytical Gradient ',Vector3D( M%y*M%z, M%x*M%z, M%x*M%y)

  call CPU_TIME(time1)
  print*, ' Total computational time was ', time1-time0,'s.'
  stop
END PROGRAM TEST_BEM
