include 'Lapl_3D.f90'
!------------------------------------------------------------------------!
!-----  (c) Mykola Yas'ko. All right reserved.                    -------!
!-----  http://nyasko.homepage.com   E-mail: nyasko@hotmail.com   -------!
!-----  This is Demo-program for solving 3D Laplace equation      -------!
!-----  by the direct boundary element method.                    -------!
!-----  Domain is interior of sphere, radius=1                    -------!
!------------------------------------------------------------------------!
integer function Laplace(N,P1,P2,P3) !! initialize object coordinates
  USE BEM_Lapl_3D
  IMPLICIT NONE
  integer, intent(in) :: N
  type(Vector3D), dimension(N), intent(in) :: P1,P2,P3
  type(Vector3D)     :: M, Grad
  integer :: i, Method=4
  real(bem) :: Pot
  real :: time0, time1

!---------- Print vector parameters ---------------------
!    PRINT *, "1 -------------"
!    WRITE(*,"(3F8.2)")  P1(1)%x, P1(1)%y, P1(1)%z
!    PRINT *, "2 -------------"
!    WRITE(*,"(3F8.2)")  P2(1)%x, P2(1)%y, P2(1)%z
!    PRINT *, "3 -------------"
!    WRITE(*,"(3F8.2)")  P3(1)%x, P3(1)%y, P3(1)%z
!    PRINT *, "4 -------------"

!---------- Allocate N boundary elements ----------------
  i=InitialiseBE(N)
!--------------------------------------------------------

!-------- Input boundary elements vertices --------------
  do i=1,N
     BE(i)%P1=P1(i)
     BE(i)%P2=P2(i)
     BE(i)%P3=P3(i)
  end do
!--------------------------------------------------------

  i = FillBE()     ! Additional computation for boundary elements

  if ( i == -7 ) then
     stop ' Bad boundary elements occured. Stop!'
  elseif ( i == -6 ) then
     stop ' Function InitialiseBE did not used. Stop!'
  end if


!-------- Boundary condition of 1nd kind  --- ----------------------
 BE%bc = 1   ! Kind of boundary condition
 BE%POT = 1._bem
!-------------------------------------------------------------------

  i = SolveBVP(Method)         ! Solve boundary-valued problem by
  Laplace = i
  return
END function Laplace

!integer function Testfun(N)
!  integer :: N
!  Testfun = N*2
!  return
!END function Testfun

  subroutine ResultF(Point, Pot, Grad)
  USE BEM_Lapl_3D
    type(Vector3D), intent(in) :: Point
    real(bem), intent(out) :: Pot
    type(Vector3D), intent(out) :: Grad
    Pot = Potential(Point)
    Grad = Gradient(Point)
  end subroutine ResultF

