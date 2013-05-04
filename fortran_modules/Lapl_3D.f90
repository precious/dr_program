module Lin_Sys      ! (c) Mykola Yas'ko http://www.uni-koblenz.de/~yasko
  implicit none     !                   E-mai:  nyasko@hotmail.com                  
  integer, parameter, public  :: bem=selected_real_kind(12,36)
  integer, parameter, public  :: dp=max(selected_real_kind(14,50),bem)
  integer, parameter, public  :: mGaussElim=0, mCoGrad=1, mBiCoGrad=2, mPCGSM=3, mPBiCGSM=4, MUnknown=200
  real(bem), parameter :: pi=3.14159265358979324_bem, qpi=1._bem/pi, qpi2=0.5_bem/pi, pid2=pi/2._bem, pi4=pi/4._bem,&
                          pi5=0.5_bem*pi, qpi4=0.25_bem/pi, pi2=pi+pi, qpi8=0.125_bem/pi, dpi=2._bem/pi
  interface GaussElim
    module procedure GAUSS_ELIM
    module procedure GAUSS_ELIM_N
  end interface GaussElim 

contains

 subroutine SolveLineSystem( A, K, X, Ierr, Method )
  integer, intent(in)   :: Method
  real(bem), INTENT(in out) :: A(:,:), K(:)
  real(bem), INTENT(in out) :: X(:)
  integer, intent(in out)      :: Ierr
   if(Method == mCoGrad) then
      call COGRAD( A, K, X, Ierr)
   elseif( Method == mBiCoGrad ) then
      call BiCoGrad( A, K, X, Ierr)
   elseif( Method == mPCGSM ) then
      call PCGSM ( A, K, X, Ierr )
   elseif( Method == mPBiCGSM ) then
      call BiPCGSM ( A, K, X, Ierr )
   else
      call GAUSS_ELIM ( A, K, X, Ierr )
   end if
   return
 end subroutine SolveLineSystem

 subroutine SolveLineSystemC( A, K, X, Ierr, Method )
  integer, intent(in)          :: Method
  real(bem), INTENT(in out)    :: A(:,:)
  complex(bem), INTENT(in out) :: K(:), X(:)
  integer, intent(in out)      :: Ierr
  real(bem), dimension(2*size(X)) :: KK, XX
  integer :: n, MaxIt
   n= size(X)
   KK = (/ real(K), aimag(K) /)
   XX = (/ real(X), aimag(X) /)
   if(Method == mCoGrad) then
      MaxIt=Ierr
      call COGRAD( A, KK(1:n), XX(1:n), Ierr)
      if( Ierr<0 ) return
      call COGRAD( A, KK(n+1:n+n), XX(n+1:n+n), MaxIt)
   elseif( Method == mBiCoGrad ) then
      MaxIt=Ierr
      call BiCoGrad( A, KK(1:n), XX(1:n), Ierr)
      if( Ierr<0 ) return
      call BiCoGrad( A, KK(n+1:n+n), XX(n+1:n+n), MaxIt)
   elseif( Method == mPCGSM ) then
      MaxIt=Ierr
      call PCGSM( A, KK(1:n), XX(1:n), Ierr)
      if( Ierr<0 ) return
      call PCGSM( A, KK(n+1:n+n), XX(n+1:n+n), MaxIt)
   elseif( Method == mPBiCGSM ) then
      MaxIt=Ierr
      call BiPCGSM( A, KK(1:n), XX(1:n), Ierr)
      if( Ierr<0 ) return
      call BiPCGSM( A, KK(n+1:n+n), XX(n+1:n+n), MaxIt)
   else
      call GAUSS_ELIM ( A, KK, XX, Ierr )
   end if
   if( Ierr>=0 ) X = CMPLX( XX(1:n), XX(n+1:n+n) )
   return
 end subroutine SolveLineSystemC

 SUBROUTINE COGRAD( A, K, X, ITR )
  IMPLICIT NONE
  integer, INTENT(in out)       :: ITR
  real(bem), INTENT(IN out)     :: A(:,:), K(:)
  real(bem), INTENT(IN OUT)     :: X(:)
  real(dp), PARAMETER           :: TINY=1E-33, EPS=1E-9
  real(dp), DIMENSION(SIZE(X))  :: P, R, ATRi, APi
  real(dp)                      :: S1, S2, ALF, BET
  integer                       :: N, IT, i, ItrMax
   N = SIZE(X)
   ItrMax=ITR
!---------------- Preconditioner ------------------------------
  do i=1,N
     ALF = MaxVal(abs(A(i,:)))
     A(i,:) = A(i,:)/ALF
     K(i) =  K(i)/ALF
  end do
!--------------------------------------------------------------
  R = K - MatMul(A,X)
  P = MatMul(R,A)
  ATRi = P
  S1 = TINY + Dot_Product(ATRi,ATRi)
  DO IT = 1,ItrMax
    APi = MatMul(A,P)
    S2 = TINY + Dot_Product(APi,APi)
    ALF = S1 / S2
    ATRi = ALF * P
    X = X + ATRi
    IF (MAXVAL(ABS(ATRi)) < EPS ) THEN
       ITR = IT
       RETURN
    END IF
    R = R - ALF * APi
    ATRi = MatMul(R,A)
    S2 = TINY + Dot_Product(ATRi,ATRi)
    BET = S2 / S1
    S1 = S2
    P = ATRi + BET * P
  END DO
  ITR = -7
  RETURN
 END SUBROUTINE COGRAD

 SUBROUTINE BiCoGrad( A, B, X, ITR )   ! Biconjugate gradient method with
  IMPLICIT NONE                       ! Jacobi preconditioner
  INTEGER, INTENT(in out)      :: ITR 
  REAL(bem), INTENT(IN)        :: A(:,:), B(:)
  REAL(bem), INTENT(IN OUT)    :: X(:)
  REAL(dp), PARAMETER          :: eps=1E-9
  REAL(dp), DIMENSION(SIZE(X)) :: r0, rt0, z0, zt0, p, pt, q, qt, div
  REAL(dp)                     :: alpha, ro, ro2, sum
  INTEGER                      :: N, IT, i, j
  N = SIZE(X)
  r0 = b - MatMul(A,x)
  rt0 = r0
!  div =  (/ (0.2*maxval(abs(A(i,:)))+abs(A(i,i)), i=1,N) /)
  div =  (/ (maxval(abs(A(i,:))), i=1,N) /)
  DO IT = 1,ITR
!----------- Jacobi preconditioner -----------------------
    z0 = r0 / div
    zt0 = rt0 / div
!---------------------------------------------------------
    ro =  Dot_Product(z0,rt0)
    if ( It == 1 ) then
       p = z0
       pt = zt0
    else
       if (abs(ro) <= tiny(1._dp) ) then
          ITR = -6  ! method fails
          return
       end if
       alpha = ro/ro2
       p = z0 + alpha * p
       pt = zt0 + alpha * pt
    end if
    q = MatMul(A,p)
    qt = MatMul(pt,A)
    alpha = ro / Dot_product(pt,q)
    z0 = alpha * p
    x = x + z0
    if ( MaxVal(abs(z0)) < eps ) then
       ITR = It
       return
     else
    end if
    r0 = r0 - alpha * q
    rt0 = rt0 - alpha * qt
    ro2 = ro
  END DO
  ITR = -7
  return
 END SUBROUTINE BiCoGrad

 SUBROUTINE PCGSM( A, b, x, Itr ) ! Preconditioned Conjugate Gradient Squared
  IMPLICIT NONE                  ! Method with Jacobi preconditioner
  integer, INTENT(in out)      :: Itr 
  REAL(bem), INTENT(IN)        :: A(:,:), B(:)
  REAL(bem), INTENT(IN OUT)    :: X(:)
  REAL(dp), PARAMETER          :: eps=1E-8, tin=1E-33
  REAL(dp), DIMENSION(SIZE(X)) :: r0,  u1, p1, q0, rt, pt, vt, div
  REAL(dp)                     :: alf, ro, ro2
  integer                      :: N, It, i, ItrMax
  N = SIZE(X)
  ItrMax=ITR
  r0 = b - MatMul(A,x)
  rt = r0
  div = (/ (maxval(abs(A(i,:))), i=1,N) /)
  do It = 1,ItrMax
    ro = Dot_Product(rt, r0)
    if(abs(ro)<tin .and. maxval(abs(r0))>eps) then
       Itr = -16
       return
    end if
    if(It==1) then
       u1 = r0
       p1 = u1
    else
       alf = ro/ro2
       u1 = r0 + alf*q0
       p1 = u1 + alf*(q0+alf*p1)
    end if
!----------- Jacobi preconditioner -----------------------
    pt = p1/div
!---------------------------------------------------------
    vt = matmul(A,pt)
    alf = ro/Dot_Product(rt,vt)
    q0 = u1 - alf*vt
!----------- Jacobi preconditioner -----------------------
    u1 = alf*(u1+q0) / div
!---------------------------------------------------------
    x = x + u1
    if ( MaxVal(abs(u1)) < eps ) then
       Itr = It
       return
     else
    end if
    r0 = r0 - matmul(A,u1)
    ro2 = ro
  end do
  ITR = -7
  return
 END SUBROUTINE PCGSM

 SUBROUTINE BiPCGSM( A, b, x, Itr ) ! BiConjugate Gradient Stabilized
  IMPLICIT NONE                    ! Method 
  integer, INTENT(in out)      :: Itr 
  REAL(bem), INTENT(IN)        :: A(:,:), B(:)
  REAL(bem), INTENT(IN OUT)    :: X(:)
  REAL(dp), PARAMETER          :: eps=1E-9, tin=1E-34
  REAL(dp), DIMENSION(SIZE(X)) :: r0, rt, p, s, t, vt
  REAL(dp)                     :: alf, bet, omeg,  ro, ro2
  integer                      :: N, It, ItrMax
  N = SIZE(X)
  ItrMax=ITR
  r0 = b - MatMul(A,x)
  rt = r0
  do It = 1,ItrMax
    ro = Dot_Product(rt, r0)
    if(abs(ro)<tin .and. maxval(abs(r0))>eps .and. It>1) then
       Itr = -16
       return
    end if
    if(It==1) then
       p = r0
    else
       bet = ro/ro2*(alf/omeg)
       p = r0 + bet*(p-omeg*vt)
    end if
    vt = matmul(A,p)
    alf = ro/Dot_Product(rt,vt)
    s = r0 - alf*vt
    if ( abs(alf)*MaxVal(abs(s)) < eps ) then
       x = x + alf*p
       Itr = It
       return
    end if
    t = matmul(A,s)
    omeg = Dot_Product(t,s)/Dot_Product(t,t)
    x = x + alf*p + omeg*s
    ro2 = ro
    r0 = s - omeg*t
  end do
  ITR = -7
  return
 END SUBROUTINE BiPCGSM

 subroutine GAUSS_ELIM ( A, Y, X, Ierror )  ! Gauss elimination method
  IMPLICIT NONE
  real(bem), INTENT(IN OUT) :: A(:,:), X(:), Y(:)
  integer, intent(out) :: Ierror
  real(bem) :: alpha, TEMP, SWAP(SIZE(X))
  integer :: i,j, n, location, location_array(1)
   n = SIZE( X )
   do i = 1, n-1
     location_array = MAXLOC( ABS(A(i:n,i)) )
     location = location_array(1)+i-1
     if (location /= i) then
       SWAP(i:n) = A(i,i:n)
       A(i,i:n) = A(location,i:n)
       A(location,i:n) = SWAP(i:n)
       TEMP = Y(i)
       Y(i) = Y(location)
       Y(location) = TEMP
     end if
     do j = i+1, n
       alpha = A(j,i) / A(i,i)
       A(j, i:n) = A( j, i:n ) - alpha * A(i,i:n)
       Y(j) = Y(j) - alpha * Y(i)
     end do
   end do
   if( any( (/ (abs(A(i,i)), i=1,n) /) <= tiny(1._bem)) ) then
        Ierror=-1
        return
   end if
   X(n) = Y(n) / A(n,n)
   do i = n-1, 1, -1
     Y(1:i) = Y(1:i) - X(i+1) * A(1:i, i+1)
     X(i) = Y(i) / A(i,i)
   end do
   Ierror=0
   return
 end subroutine GAUSS_ELIM

 subroutine GAUSS_ELIM_N ( A, Y, X, Ierror )
  IMPLICIT NONE
  real(bem), INTENT(IN OUT) :: A(:,:), X(:,:), Y(:,:)
  integer, intent(out) :: Ierror
  real(bem) :: alpha, TEMP(SIZE(X,2)), SWAP(SIZE(X,1))
  integer :: i,j,k,n, location, location_array(1)
   n = SIZE( X, 1 )
   k = SIZE( X, 2 )
   do i = 1, n-1
     location_array = MAXLOC( ABS(A(i:n,i)) )
     location = location_array(1)+i-1
     if (location /= i) then
       SWAP(i:n) = A(i,i:n)
       A(i,i:n) = A(location,i:n)
       A(location,i:n) = SWAP(i:n)
       TEMP = Y(i,:)
       Y(i,:) = Y(location,:)
       Y(location,:) = TEMP
     end if
     do j = i+1, n
       alpha = A(j,i) / A(i,i)
       A(j, i:n) = A( j, i:n ) - alpha * A(i,i:n)
       Y(j,:) = Y(j,:) - alpha * Y(i,:)
     end do
   end do
   if( any( (/ (abs(A(i,i)), i=1,n) /) <= tiny(1._bem)) ) then
        Ierror=-1
        return
   end if
   X(n,:) = Y(n,:) / A(n,n)
   do i = n-1, 1, -1
     ForAll(j=1:k) Y(1:i,j) = Y(1:i,j) - X(i+1,j) * A(1:i, i+1)
     X(i,:) = Y(i,:) / A(i,i)
   end do
   Ierror=0
   return
 end subroutine GAUSS_ELIM_N

 function Determinant( A ) result(Pr)
  IMPLICIT NONE
  real(bem), INTENT(IN OUT) :: A(:,:)
  real(bem)                 :: Pr
  real(bem) :: SWAP(UBOUND( A, 1 )), alpha
  integer :: i,j, n
   n = UBOUND( A, 1 )
   if(abs(A(1,1)) < 1e-9_bem ) A(1,1)= sign(1e-9_bem, A(1,1))
   Pr = A(1,1)
   do i = 1, n-1
     if(abs(A(i,i)) < 1e-9_bem ) A(i,i)= sign(1e-9_bem,A(i,i))
     do j = i+1, n
       alpha = A(j,i) / A(i,i)
       A(j, i+1:n) = A( j, i+1:n) - alpha * A(i,i+1:n)
     end do
     Pr = Pr * A( i+1, i+1 )
     if( abs(Pr) < tiny(0._bem) ) return
   end do
   return
 end function Determinant

 function Det2( A )  result(Pr)
  IMPLICIT NONE
  real(bem), INTENT(IN) :: A(:,:)
  real(bem)             :: Pr
  real(bem) :: B(UBOUND(A,1),UBOUND(A,2))
  integer :: i,j, n
   n = UBOUND( A, 1 )
   B(1,:) = A(1,:)
   B(2:n,1) = A(2:n,1)/B(1,1)
   Pr = B(1,1)
   do i = 2, n
     do j = 2, i-1
       B(i,j) = (A(i,j)-Dot_Product(B(i,1:j-1),B(1:j-1,j)))/B(j,j)
     end do
     do j=i,n
       B(i,j) = A(i,j)-Dot_Product(B(i,1:i-1),B(1:i-1,j))
     end do
     Pr = Pr*B(i,i)
     if( abs(Pr) < tiny(0._bem) ) return
   end do
   return
 end function Det2

 function Det2m(A) result(Pr)
  IMPLICIT NONE
  real(bem), intent(in out) :: A(:,:)
  real(bem)             :: Pr
  integer :: i,j, n
   n = UBOUND(A,1)
   A(2:n,1) = A(2:n,1)/A(1,1)
   Pr = A(1,1)
   do i = 2, n
     do j = 2, i-1
       A(i,j) = (A(i,j)-Dot_Product(A(i,1:j-1),A(1:j-1,j)))/A(j,j)
     end do
     do j=i,n
       A(i,j) = A(i,j)-Dot_Product(A(i,1:i-1),A(1:i-1,j))
     end do
     Pr = Pr*A(i,i)
     if( abs(Pr) < tiny(0._bem) ) return
   end do
   return
 end function Det2m
end module Lin_Sys

module BEM_3D
 use Lin_Sys
   IMPLICIT NONE

   TYPE, public  :: Vector3D
     real(bem)   ::  x,y,z
   END TYPE Vector3D

   interface operator (+)
     module procedure Plus00
   end interface 

   interface operator (-)
     module procedure Minus0, Minus00
   end interface

   interface operator (*)
     module procedure Mult0c, Multc0
   end interface 

   interface operator (/)
     module procedure Div0c0
   end interface 

   interface operator (**)
     module procedure Deg0
   end interface 

   interface operator (.sp.)
     module procedure ScalProd00
   end interface

   interface operator (.vp.)
     module procedure VecProd00
   end interface

   interface Length
     module procedure Length0
   end interface 

   private Plus00, Minus0, Minus00, Mult0c, Multc0, Div0c0, Deg0, ScalProd00, VecProd00, Length0

CONTAINS

   elemental function Plus00(V1,V2)
    type(Vector3D), intent(in) :: V1,V2
    type(Vector3D)             :: Plus00
     Plus00 = Vector3D( V1%x+V2%x, V1%y+V2%y, V1%z+V2%z )
     return
   end function Plus00

   elemental function Minus0(V1)
    type(Vector3D), intent(in) :: V1
    type(Vector3D)             :: Minus0
     Minus0 = Vector3D( -V1%x, -V1%y, -V1%z )
     return
   end function Minus0
   elemental function Minus00(V1,V2)
    type(Vector3D), intent(in) :: V1,V2
    type(Vector3D)             :: Minus00
     Minus00 = Vector3D( V1%x-V2%x, V1%y-V2%y, V1%z-V2%z )
     return
   end function Minus00

   elemental function Mult0c(V1,V2)
    type(Vector3D), intent(in) :: V1
    real(bem), intent(in)      :: V2
    type(Vector3D)             :: Mult0c
     Mult0c = Vector3D(V1%x*V2,V1%y*V2,V1%z*V2)
     return
   end function Mult0c
   elemental function Multc0(V1,V2)
    real(bem), intent(in)      :: V1
    type(Vector3D), intent(in) :: V2
    type(Vector3D)             :: Multc0
     Multc0 = Vector3D(V1*V2%x,V1*V2%y,V1*V2%z)
     return
   end function Multc0

   elemental function Div0c0(V1,V2)
    type(Vector3D), intent(in) :: V1
    real(bem), intent(in)      :: V2
    type(Vector3D)             :: Div0c0
     Div0c0 = Vector3D(V1%x/V2,V1%y/V2,V1%z/V2)
     return
   end function Div0c0

   elemental function Deg0(V,n)
    type(Vector3D), intent(in) :: V
    integer, intent(in) :: n
    real(bem)           :: Deg0
     Deg0 = V%x**n+V%y**n+V%z**n
     return
   end function Deg0

   elemental function Length0(V)
    type(Vector3D), intent(in) :: V
    real(bem)                 :: Length0
     Length0 = sqrt(V%x**2+V%y**2+V%z**2)
     return
   end function Length0

   elemental function ScalProd00(V1,V2)
    type(Vector3D), intent(in) :: V1,V2
    real(bem)                 :: ScalProd00
     ScalProd00 = V1%x*V2%x+V1%y*V2%y+V1%z*V2%z
     return
   end function ScalProd00

   elemental function VecProd00(V1,V2)
    type(Vector3D), intent(in) :: V1,V2
    type(Vector3D)             :: VecProd00
     VecProd00 = Vector3D(V1%y*V2%z-V1%z*V2%y, V1%z*V2%x-V1%x*V2%z, V1%x*V2%y-V1%y*V2%x )
     return
   end function VecProd00
end module BEM_3D

MODULE BEM_Lapl_3D
   USE BEM_3D
   IMPLICIT NONE
   TYPE, public :: BoundaryElement
      type(Vector3D)   :: P1,P2,P3,P, s1,s2,s3, V, Vt, n
      real(bem)        :: L12,L13,L23,SQ, POT,Vn,g
      integer          :: bc
   END TYPE BoundaryElement

   real(bem), allocatable :: AA(:,:), BB(:), CC(:)
   TYPE(BoundaryElement),  pointer, save :: BE(:)
   TYPE(BoundaryElement), allocatable, target :: BEO(:)
   integer, public    :: NBE=0, NEQ=0
   integer, private   :: i=908996950, j, State=0, k_sym=0, kx_sym=0, ky_sym=0, kz_sym=0
   real(bem), parameter, private :: EPS=1E-5
   real(bem), private :: x_sym, y_sym, z_sym
   type(Vector3D), private :: VsAdd=Vector3D(0.0_bem,0.0_bem,0.0_bem)
   logical, private        :: VS_Add=.false.

CONTAINS

   integer function InitialiseBE(N_BE, N_EQ)
     integer, intent(in)            :: N_BE
     integer, intent(in), optional  :: N_EQ
      IF ( present(N_EQ) ) then
          NEQ = N_EQ
      ELSE
          NEQ = N_BE
      END IF
      NBE = N_BE
      IF ( ALLOCATED( BEO ) ) DEALLOCATE( BEO )
      ALLOCATE ( BEO(NBE), STAT=I )
      if ( I == 0) then
        BE => BEO
        State = 1
        InitialiseBE = NBE
      else
        State = -1
        InitialiseBE = 0
      end if
      return
   END function InitialiseBE

   integer function FillBE()
      if( State < 1 ) then
         FillBE = -1
         return
      end if
      BE%P = BE%P1+BE%P2+BE%P3
      BE%P = BE%P/3._bem
      BE%s1 = BE%P2-BE%P1
      BE%L12 = Length(BE%s1)
      BE%s1 = BE%s1 / BE%L12
      BE%s2 = BE%P3-BE%P2
      BE%L23 = Length(BE%s2)
      BE%s2 = BE%s2 / BE%L23
      BE%s3 = BE%P1-BE%P3
      BE%L13 = Length(BE%s3)
      BE%s3 = BE%s3 / BE%L13
      BE%n = (BE%P2-BE%P1) .vp. (BE%P3-BE%P1)
      BE%SQ = Length( BE%n )
      BE%n =  BE%n / BE%SQ
      BE%SQ = 0.5_bem*BE%SQ
      BE%Vt = Vector3D(0.0_bem,0.0_bem,0.0_bem)
      State = 2
      FillBE = size(BE)
      return
   END Function FillBE

   integer function SetSymmetry(kx, x, ky, y, kz, z)
    integer, intent(in) :: kx,ky,kz ! Kind of symmetry: 1 - (Vn=0),  -1 - (Pot=0), 0 - no symmetry
    real(bem), intent(in) :: x,y,z  ! Planes of symmetry
     if(max(abs(kx), abs(ky), abs(kz))>1) then
        SetSymmetry = -21 ! Bad kind of symmetry 
        return
     end if
     k_sym = 0
     kx_sym=kx
     if(kx/=0) k_sym=k_sym+1
     ky_sym=ky
     if(ky/=0) k_sym=k_sym+1
     kz_sym=kz
     if(kz/=0) k_sym=k_sym+1
     x_sym=x
     y_sym=y
     z_sym=z
     SetSymmetry = k_sym
   end function SetSymmetry
   
  subroutine LineCoef( Point, GdFn, GFi )
   implicit none
   type(Vector3D), intent(in)           :: Point
   real(bem), dimension(:), intent(out) :: GdFn, GFi
   type(Vector3D) :: Ps(4)
   integer        :: i, N, Ns(4)
   real(bem), dimension(size(GFi))      :: GdFn1, GFi1
    call OneLineCoef( Point, GdFn, GFi )
    if(k_sym==0) return
    N = 1
    Ns(N)=1
    Ps(N) = Point
    if(kx_sym /= 0) then
       N = 2
       Ps(N) = Vector3D( 2*x_sym-Point%x, Point%y, Point%z )
       Ns(N)=kx_sym
       call OneLineCoef( Ps(N), GdFn1, GFi1 )
       GFi = GFi + kx_sym*GFi1
       GdFn = GdFn + kx_sym*GdFn1
    end if
    if(ky_sym /= 0) then
      do i=1,N
       Ps(N+i) = Vector3D( Ps(i)%x, 2*y_sym-Ps(i)%y, Ps(i)%z )
       Ns(N+i) = (Ns(i)*ky_sym)
       call OneLineCoef( Ps(N+i), GdFn1, GFi1 )
       GFi = GFi + Ns(N+i)*GFi1
       GdFn = GdFn + Ns(N+i)*GdFn1
      end do
      N = N+N
    end if
    if(kz_sym /= 0) then
      do i=1,N
       call OneLineCoef(Vector3D(Ps(i)%x,Ps(i)%y,2*z_sym-Ps(i)%z), GdFn1, GFi1 )
       GFi = GFi + (Ns(i)*kz_sym)*GFi1
       GdFn = GdFn + (Ns(i)*kz_sym)*GdFn1
      end do
    end if
  contains
   SUBROUTINE OneLineCoef( Point, GdFn, GFi )
     type(Vector3D), INTENT(IN)            :: Point
     real(bem), dimension(:), intent(out) :: GdFn, GFi
     type(Vector3D) :: R, R1,R2,R3, R1R2, R1R3, R2R3, Ra, Rb, P
     real(bem)      :: RR, sq_RR, Rn, Rn2, a1,a2,a3, b1,b2, s1,s2,s3, &
                       p1,p2,p3, q1,q2,q3
     integer        :: j
     do j=1,NBE
        P = Point
        R = BE(j)%P-P
        Rn = R .sp. BE(j)%n
        if ( abs(Rn) < Eps ) then
           P = Point - BE(j)%n*(Eps-Rn)
           R = BE(j)%P-P
           Rn = Eps
        end if
        Rn2 = Rn*Rn
        RR = R**2
        sq_RR = sqrt(RR)
        if ( RR >  121._bem*BE(j)%SQ ) then        ! Linear       1 node
           GdFn(j) = qpi4 / sq_RR * BE(j)%SQ
           GFi(j)  = qpi4 * Rn /(sq_RR*RR) * BE(j)%SQ
        else if ( RR > 64._bem*BE(j)%SQ ) then        !     3 nodes
           R = (0.5_bem*(BE(j)%P1+BE(j)%P2)) - Point
           a1 = 1./Length(R)
           s1 = a1**3
           b1 = a1
           b2 = s1
           R = (0.5_bem*(BE(j)%P2+BE(j)%P3)) - Point
           a1 = 1./Length(R)
           s1 = a1**3
           b1 = b1 + a1
           b2 = b2 + s1
           R = (0.5_bem*(BE(j)%P3+BE(j)%P1)) - Point
           a1 = 1/Length(R)
           s1 = a1**3
           b1 = b1 + a1
           b2 = b2 + s1
           GdFn(j) = qpi4 * b1 * BE(j)%SQ / 3.0_bem
           GFi(j)  = qpi4 * Rn * b2 * BE(j)%SQ / 3.0_bem
        else
           R1 = BE(j)%P1 - P
           R2 = BE(j)%P2 - P
           R3 = BE(j)%P3 - P
           R1R2 = R1 .vp. R2
           R1R3 = R1 .vp. R3
           R2R3 = R2 .vp. R3
           Ra = R1R2 .vp. R1
           Rb = R1R3 .vp. R1
           q1 = (Ra .sp. Rb)/sqrt((Ra.sp.Ra)*(Rb.sp.Rb))
           q1 = acos( max(-1.0_bem, min(1.0_bem,q1)) )
           Ra = R2R3 .vp. R2
           Rb = -R1R2 .vp. R2
           q2 = (Ra .sp. Rb)/sqrt((Ra.sp.Ra)*(Rb.sp.Rb))
           q2 = acos( max(-1.0_bem, min(1.0_bem,q2)) )
           Ra = -R1R3 .vp. R3
           Rb = -R2R3 .vp. R3
           q3 = (Ra .sp. Rb)/sqrt((Ra.sp.Ra)*(Rb.sp.Rb))
           q3 = acos( max(-1.0_bem, min(1.0_bem,q3)) ) 
           GFi(j) = qpi4 * sign( q1+q2+q3 - pi , Rn)
           s1=sqrt(R1 .sp. R1)
           s2=sqrt(R2 .sp. R2)
           s3=sqrt(R3 .sp. R3)
           q1 = R1 .sp. BE(j)%s1
           q2 = R2 .sp. BE(j)%s2
           q3 = R3 .sp. BE(j)%s3
           Ra=(log((BE(j)%L12+q1+s2)/(q1+s1))*(R1.vp.BE(j)%s1)) &
             +(log((BE(j)%L23+q2+s3)/(q2+s2))*(R2.vp.BE(j)%s2)) &
             +(log((BE(j)%L13+q3+s1)/(q3+s3))*(R3.vp.BE(j)%s3))
           GdFn(j)= qpi4*(Ra.sp.BE(j)%n) - Rn*GFi(j)
        end if
     end do
     RETURN
   END SUBROUTINE OneLineCoef
  END SUBROUTINE LineCoef

   function SolveBVP( Method ) result (Error)
     integer, intent(in) :: Method
     integer             :: Error
     integer        :: i,j
     real(bem)      :: GdFn(NBE), GFi(NBE)
     type(Vector3D) :: VFi(NBE), VdFn(NBE), Vn
   VS_Add=.false.
   if( State /= 2 ) then
     Error = -5   ! Function FillBE did not used
     return
   elseif( .not. associated(BE) ) then
     Error = -4   ! Boundary Elements are out of memory
     return
   elseif(any(BE%bc<1) .or. any(BE%bc>5) ) then
     Error = -10  ! Unknown boundary condition
     return
   elseif( NBE<3) then
     Error = -2   ! Small number of boundary elements
     return
   end if
   ALLOCATE ( AA(NEQ,NEQ), BB(NEQ), CC(NEQ), stat=i)
   if(i /= 0) Then
      Error=-6
      return
   end if
     do i=1,NBE
       if(BE(i)%bc==5) then
          VS_Add=.true.
          VsAdd=Vector3D(0.0_bem,0.0_bem,0.0_bem)
          call GradLine( BE(i)%P, VFi, VdFn )
          BB(i) = -(VsAdd.sp.BE(i)%n)
          VS_Add=.false.
          do j=1,NBE
            if(BE(j)%bc==1.or.BE(j)%bc==5) then
              AA(i,j) = VdFn(j).sp.BE(i)%n
              BB(i) = BB(i) - (VFI(j).sp.BE(i)%n)*BE(j)%POT
            elseif(BE(j)%bc==2) then
              AA(i,j) = VFI(j).sp.BE(i)%n
              BB(i) = BB(i) - (VdFn(j).sp.BE(i)%n)*BE(j)%Vn
            elseif(BE(j)%bc==3) then
              AA(i,j) = (VFI(j) + (BE(j)%Vn*VdFn(j))).sp.BE(i)%n
              BB(i) = BB(i) - (VdFn(j).sp.BE(i)%n)*BE(j)%g
            end if
          end do
          AA(i,i) = AA(i,i) - 1._bem
       else
          BB(i) = 0.0_bem
          call LineCoef( BE(i)%P, GdFn, GFi )
          GFi(i) = GFi(i) - 1._bem
          do j=1,NBE
            if(BE(j)%bc==1.or.BE(j)%bc>=5) then
              AA(i,j) = GdFn(j)
              BB(i) = BB(i) - GFI(j)*BE(j)%POT
            elseif(BE(j)%bc==2) then
              AA(i,j) = GFI(j)
              BB(i) = BB(i) - GdFn(j)*BE(j)%Vn
            elseif(BE(j)%bc==3) then
              AA(i,j) =  GFI(j) + BE(j)%Vn*GdFn(j)
              BB(i) = BB(i) - GdFn(j)*BE(j)%g
            end if
          end do
       end if
     end do
     I = NBE
     CALL SolveLineSystem( AA, BB, CC, I,  MOD(Method,10) ) ! Solve system of linear equation
     Error = I
     if( Error >= 0 ) then
      do i=1,NBE
        if(BE(i)%bc==1.or.BE(i)%bc==5) then
          BE(i)%Vn = CC(i)
        else if(BE(i)%bc == 2) then
          BE(i)%POT = CC(i)
        else if(BE(i)%bc == 3) then
          BE(i)%POT = CC(i)
          BE(i)%Vn = BE(i)%Vn*CC(i) + BE(i)%g
        end if
      end do 
     end if
     DEALLOCATE (AA, BB, CC)
     if(Method<10) then
      do i=1,NBE
        if(BE(i)%bc<5) then
          Vn = BE(i)%Vn*BE(i)%n
          BE(i)%Vt = (2._bem*(Gradient(BE(i)%P)-Vn))
          BE(i)%V = BE(i)%Vt + Vn
        else
          BE(i)%V = BE(i)%Vt + (BE(i)%Vn*BE(i)%n)
        end if
      end do 
      VS_Add=.true.
     end if
     RETURN
   END function SolveBVP

   FUNCTION Potential( Point )
     real(bem)                 :: Potential
     type(Vector3D), INTENT(IN) :: Point
     real(bem), DIMENSION(NBE) :: GdFn, GFi
       CALL LineCoef( Point, GdFn, GFi )
       Potential = dot_product(GdFn,BE%Vn) + dot_product(GFi,BE%POT)
      RETURN
   END FUNCTION Potential

   FUNCTION Gradient( Point )
     type(Vector3D)             :: Gradient
     type(Vector3D), INTENT(IN) :: Point
     type(Vector3D), DIMENSION(NBE) :: GFi, GdFn
       VsAdd=Vector3D(0.0_bem,0.0_bem,0.0_bem)
       call GradLine( Point, GFi, GdFn )
       Gradient = VSAdd+Vector3D(dot_product(GdFn%x,BE%Vn)+dot_product(GFi%x,BE%Pot), &
                           dot_product(GdFn%y,BE%Vn)+dot_product(GFi%y,BE%Pot), &
                           dot_product(GdFn%z,BE%Vn)+dot_product(GFi%z,BE%Pot) )
      RETURN
   END FUNCTION Gradient

  subroutine GradLine( Point, GFi, GdFn )
   type(Vector3D), INTENT(IN) :: Point
   type(Vector3D), dimension(Nbe),intent(out) :: GFi, GdFn
   type(Vector3D), dimension(Nbe) :: GFi1, GdFn1
   type(Vector3D) :: Ps(4)
   integer        :: i, N, Ns(4)
    call  OneGradLine(Point, GFi, GdFn)
    if(k_sym==0) return
    N = 1
    Ps(N) = Point
    Ns(1) = 1
    if(kx_sym /= 0) then
       N = 2
       Ps(N) = Vector3D( 2*x_sym-Point%x, Point%y, Point%z )
       Ns(N) = kx_sym
       call OneGradLine(Ps(N), GFi1, GdFn1)
       GFi1%x = -GFi1%x
       GdFn1%x = -GdFn1%x
       if(kx_sym==1) then
         GFi = GFi + GFi1
         GdFn = GdFn + GdFn1
       else
         GFi = GFi - GFi1
         GdFn = GdFn - GdFn1
       end if
    end if
    if(ky_sym /= 0) then
      do i=1,N
       Ps(N+i) = Vector3D( Ps(i)%x, 2*y_sym-Ps(i)%y, Ps(i)%z )
       Ns(N+i) = Ns(i)*ky_sym
       call OneGradLine(Ps(N+i), GFi1, GdFn1)
       if(kx_sym/=0 .and. i==2) then
         GFi1%x = -GFi1%x
         GdFn1%x = -GdFn1%x
       end if
       GFi1%y = -GFi1%y
       GdFn1%y = -GdFn1%y
       if(Ns(N+i) == 1) then
         GFi = GFi + GFi1
         GdFn = GdFn + GdFn1
       else
         GFi = GFi - GFi1
         GdFn = GdFn - GdFn1
       end if
      end do
      N = N+N
    end if
    if(kz_sym /= 0) then
      do i=1,N
       call OneGradLine(Vector3D(Ps(i)%x,Ps(i)%y,2*z_sym-Ps(i)%z), GFi1, GdFn1)
       if(N==2 .and. i==2) then
        if(kx_sym/=0) then
         GFi1%x = -GFi1%x
         GdFn1%x = -GdFn1%x
        else if(ky_sym/=0) then
         GFi1%y = -GFi1%y
         GdFn1%y = -GdFn1%y
        end if
       end if
       if(N==4) then
         if(i==2 .or. i==4) then
           GFi1%x = -GFi1%x
           GdFn1%x = -GdFn1%x
         end if
         if(i>=3) then
           GFi1%y = -GFi1%y
           GdFn1%y = -GdFn1%y
         end if
       end if
       GFi1%z = -GFi1%z
       GdFn1%z = -GdFn1%z
       if(Ns(i)*kz_sym==1) then
         GFi = GFi + GFi1
         GdFn = GdFn + GdFn1
       else
         GFi = GFi - GFi1
         GdFn = GdFn - GdFn1
       end if      
      end do
    end if
  return
  contains
   subroutine OneGradLine(Point, GFi, GdFn)
     type(Vector3D), INTENT(IN) :: Point
     type(Vector3D), dimension(Nbe), intent(out) :: GFi, GdFn
     type(Vector3D) :: R, R1,R2,R3, R1R2, R1R3, R2R3, Ra, Rb, Rd, P, Rs1, Rs2, Rs3
     real(bem)      :: RR, sq_RR, Rn, Tet, a1,a2,a3, s1,s2,s3, q1,q2,q3, l1,l2,l3
     integer        :: j
     do j=1,NBE
        P = Point
        R = BE(j)%P-P
        Rn = R .sp. BE(j)%n
        if ( abs(Rn) < Eps ) then
           P = Point - BE(j)%n*(Eps-Rn)
           R = BE(j)%P-P
           Rn = Eps
        end if
        RR = R**2
        sq_RR = sqrt(RR)
        if ( RR >   64._bem*BE(j)%SQ ) then        ! Linear       1 node
           a1 = qpi4*BE(j)%SQ/(RR*sq_RR)
           GdFn(j) =  a1 * R
           GFi(j) =   (-a1)*(BE(j)%n-((3.0_bem*Rn/RR)*R))
        else                 ! Analytical formulae
           R1 = BE(j)%P1 - P
           R2 = BE(j)%P2 - P
           R3 = BE(j)%P3 - P
           R1R2 = R1 .vp. R2
           R1R3 = R1 .vp. R3
           R2R3 = R2 .vp. R3
           Ra = R1R2 .vp. R1
           Rb = R1R3 .vp. R1
           a1 = (Ra .sp. Rb)/(Length(Ra)*Length(Rb))
           a1 = acos( max(-1.0_bem, min(1.0_bem,a1)) )
           Ra = R2R3 .vp. R2
           Rb = -R1R2 .vp. R2
           a2 = (Ra .sp. Rb)/(Length(Ra)*Length(Rb))
           a2 = acos( max(-1.0_bem, min(1.0_bem,a2)) )
           Ra = -R1R3 .vp. R3
           Rb = -R2R3 .vp. R3
           a3 = (Ra .sp. Rb)/(Length(Ra)*Length(Rb))
           a3 = acos( max(-1.0_bem, min(1.0_bem,a3)) ) 
           Tet = qpi4 * sign( a1+a2+a3 - pi , Rn)  ! Gap here
           a1 = R1 .sp. R1
           a2 = R2 .sp. R2
           a3 = R3 .sp. R3
           l1=sqrt(a1)
           l2=sqrt(a2)
           l3=sqrt(a3)
           Rs1 = R1 .vp. BE(j)%s1
           Rs2 = R2 .vp. BE(j)%s2
           Rs3 = R3 .vp. BE(j)%s3
           q1 = R1 .sp. BE(j)%s1
           q2 = R2 .sp. BE(j)%s2
           q3 = R3 .sp. BE(j)%s3
           s1 = ((BE(j)%L12+q1)/l2-q1/l1)/(a1-q1*q1)
           s2 = ((BE(j)%L23+q2)/l3-q2/l2)/(a2-q2*q2)
           s3 = ((BE(j)%L13+q3)/l1-q3/l3)/(a3-q3*q3)
           GFi(j) =  qpi4*( (s1*Rs1) + (s2*Rs2) + (s3*Rs3) )
           Rd=(log((BE(j)%L12+q1+l2)/(q1+l1))*BE(j)%s1)+(log((BE(j)%L23+q2+l3)/(q2+l2))*BE(j)%s2) &
               +(log((BE(j)%L13+q3+l1)/(q3+l3))*BE(j)%s3)
           GdFn(j) =   (Tet*BE(j)%n)-(qpi4*(Rd.vp.BE(j)%n))
           if(VS_Add .and. RR <26.0_bem*BE(j)%SQ) then
             R1 = BE(j)%P1 - BE(j)%P
             R2 = BE(j)%P2 - BE(j)%P
             R3 = BE(j)%P3 - BE(j)%P
             Ra=BE(j)%Vt
             VsAdd=VsAdd+(qpi4*( ((s1*(Ra.sp.R1)+(l1-(a1+q1*BE(j)%L12)/l2)/(a1-q1*q1)*(Ra.sp.BE(j)%s1))*Rs1) &
                               + ((s2*(Ra.sp.R2)+(l2-(a2+q2*BE(j)%L23)/l3)/(a2-q2*q2)*(Ra.sp.BE(j)%s2))*Rs2) &
                               + ((s3*(Ra.sp.R3)+(l3-(a3+q3*BE(j)%L13)/l1)/(a3-q3*q3)*(Ra.sp.BE(j)%s3))*Rs3))) &
                               -((Ra.sp.GdFn(j))*BE(j)%n) + (Tet*Ra)
           end if
        end if
     end do
     RETURN
   END subroutine OneGradLine
  END subroutine GradLine

   subroutine Result( Point, Pot, Grad )
     type(Vector3D), INTENT(IN) :: Point
     real(bem), intent(out)     :: Pot
     type(Vector3D), intent(out) :: Grad
       Pot = Potential(Point)
       Grad = Gradient(Point)
       return
   end subroutine Result

   integer function InputBND3(FileName,ascend)
    character(*), intent(in) :: FileName
    logical, intent(in), optional :: ascend  ! if .true. the exterior problem
    character(4) :: SIGN
    character(400) :: str
    integer :: N, i, j, l, kx, ky, kz
    logical :: exter=.false.
    real(bem) :: x,y,z
     kx=0
     ky=0
     kz=0
     open(17,file=FileName, status='OLD', position='REWIND', iostat=i)
     if(i/=0) then
       InputBND3=-25  ! Bad FileName
       return
     end if
     read(17,'(a)',iostat=i) str
     if(i==0) then
       read(str,*,iostat=i) SIGN, N
       j=index(str,'x=') !    x
       if(j>4) then
          read(str(j+2:),*,iostat=l) x
          if(l/=0) then
            InputBND3=-219  ! Bad symmetry plane
            return
          end if
          kx=-1
       else
          j=index(str,'X=')
          if(j>4) then
            read(str(j+2:),*,iostat=l) x
            if(l/=0) then
              InputBND3=-219  ! Bad symmetry plane
             return
            end if
            kx=1
          end if
       end if
       j=index(str,'y=') !    y
       if(j>4) then
          read(str(j+2:),*,iostat=l) y
          if(l/=0) then
            InputBND3=-219  ! Bad symmetry plane
            return
          end if
          ky=-1
       else
          j=index(str,'Y=')
          if(j>4) then
            read(str(j+2:),*,iostat=l) y
            if(l/=0) then
              InputBND3=-219  ! Bad symmetry plane
              return
            end if
            ky=1
          end if
       end if
       j=index(str,'z=') !    z
       if(j>4) then
          read(str(j+2:),*,iostat=l) z
          if(l/=0) then
            InputBND3=-219  ! Bad symmetry plane
            return
          end if
          kz=-1
       else
          j=index(str,'Z=')
          if(j>4) then
            read(str(j+2:),*,iostat=l) z
            if(l/=0) then
              InputBND3=-219  ! Bad symmetry plane
              return
            end if
            kz=1
          end if
       end if
     end if
     if(i/=0) then
       InputBND3=-216  ! Bad file
       return
     elseif(SIGN/='BND3') then
       InputBND3=-217  ! Bad file signature
       return
     elseif(N<4) then
       InputBND3=-218  ! Bad number of boundary elements
       return
     end if
     i = InitialiseBE(N)
     if(i/=N) then
       InputBND3=i
       return
     end if
     exter=.false.
     if(present(ascend)) then
       if(ascend) exter=.true.
     end if
     do j=1,N
       read(17,'(A400)', iostat=i) str
       if(i==0) then
       if(exter) then
          read(str,*, iostat=i) BE(j)%P1, BE(j)%P3, BE(j)%P2, BE(j)%bc
          if(i/=0)  read(str,*, iostat=i) BE(j)%P1, BE(j)%P3, BE(j)%P2
       else
          read(str,*, iostat=i) BE(j)%P1, BE(j)%P2, BE(j)%P3, BE(j)%bc
          if(i/=0)   read(str,*, iostat=i) BE(j)%P1, BE(j)%P2, BE(j)%P3
       end if
       end if
       if(i/=0) then
         InputBND3=-1000001-j  ! Bad j line
         return
       end if       
     end do
     close(17)
     if(abs(kx)+abs(ky)+abs(kz)/=0) i=SetSymmetry(kx,x, ky,y, kz,z)
     InputBND3=N
     return
   end function InputBND3

END MODULE BEM_Lapl_3D
