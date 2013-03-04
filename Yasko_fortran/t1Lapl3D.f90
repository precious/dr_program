include 'Lapl_3D.f90'
!------------------------------------------------------------------------!
!-----  (c) Mykola Yas'ko. All right reserved.                    -------!
!-----  http://nyasko.homepage.com   E-mail: nyasko@hotmail.com   -------!
!-----  This is Demo-program for solving 3D Laplace equation      -------!
!-----  by the direct boundary element method.                    -------!
!-----  Domain is interior of sphere, radius=1                    -------!
!------------------------------------------------------------------------!
PROGRAM TEST_BEM       ! Sphere, 512 BE, exterior problem
  USE BEM_Lapl_3D      ! Analytical solution  \phi = x
  IMPLICIT NONE
  integer, parameter :: N = 2048 ! Number of boundary elements. 
  type(Vector3D)     :: Point, M=Vector3D(0.0_bem,0.0_bem,0.0_bem), Grad
  integer :: i
  real(bem) :: Pot
  real :: time0, time1

  call CPU_TIME(time0)      ! Remember start time

!---------- Allocate N boundary elements ---------
  i=InitialiseBE(N)
  if( i /= N ) stop ' Out of memory!'
!-------------------------------------------------

!-------- Input boundary elements vertices -------------
  open(7,file="s2048.dat", position="rewind", status='old', iostat=i)
  if(i/=0) stop ' Bad or missing data file: s2048.dat'
  do i=1,N
     read(7,*) BE(i)%P1, BE(i)%P3, BE(i)%P2 ! Exterior problem in this order
  end do
  close(7)
!-------------------------------------------------------

  i = FillBE()     ! Additional computation for boundary elements

  if ( i == -7 ) then
     stop ' Bad boundary elements occured. Stop!'
  elseif ( i == -6 ) then
     stop ' Function InitialiseBE did not used. Stop!'
  end if

!-------- Boundary condition of 2nd kind  --- ----------------------
  BE%bc = 2   ! Kind of boundary condition
  BE%Vn =  -((BE%P-M).sp.BE%n)/Length(BE%P-M)**3    ! Normal derivative; Analytical solution U=x
!-------------------------------------------------------------------

!-------- Boundary condition of 1nd kind  --- ----------------------
! BE%bc = 1  ! Kind of boundary condition
! BE%POT =  1._bem/Length(BE%P-M)   ! Analytical solution
!-------------------------------------------------------------------

  print*, 'Exterior boundary-valued problem with', size(BE),' boundary elements'

  i = SolveBVP(mPBiCGSM)         ! Solve boundary-valued problem by 
                                 ! conjugate gradient method 
  if (i < 0 ) then             
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
  print*, '   Max error =', MaxVal(abs(BE%POT-1._bem/Length(BE%P-M))), &
                            MaxVal(abs(BE%Vn+((BE%P-M).sp.BE%n)/Length(BE%P-M)**3))

  Point = Vector3D( 1.5_bem, 1.4_bem, 1.3_bem)
!  call Result( Point, Pot, Grad)
  Pot = Potential(Point)
  print *, ' Point ', Point
  print *, ' Numerical Potential ', Pot
  print *, ' AnalyticalPotential ', 1._bem/Length(Point-M)
  print *, ' Numerical  Gradient ', Gradient(Point)
  print *, ' Analytical Gradient ', -(Point-M)/Length(Point-M)**3 

  call CPU_TIME(time1)
  print*, ' Total computational time was ', time1-time0,'s.'
  stop
END PROGRAM TEST_BEM
