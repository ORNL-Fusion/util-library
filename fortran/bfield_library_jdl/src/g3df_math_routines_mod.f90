!-----------------------------------------------------------------------------
!+ Contains special math functions for g3df
!-----------------------------------------------------------------------------
Module g3df_math_routines_mod
!
! Description:
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     04/19/2011  Ported from PENTA.  JL
! Author(s): J. Lore 07/2009 - 4/19/2011

Implicit None

Contains

!!$!-----------------------------------------------------------------------------
!!$!+ Fits a polynomial to data
!!$!-----------------------------------------------------------------------------
!!$Function polyn_fit(x,y,n,degree) & 
!!$Result(p)
!!$!
!!$! JDL 5/24/2013
!!$Use kind_mod, Only: int32, real64
!!$Implicit None
!!$Real(real64), Intent(in), Dimension(n) :: x,y
!!$Integer(int32), Intent(in) :: n,degree
!!$Real(real64), Dimension(degree+1) :: p
!!$
!!$Real(real64), Dimension(n,degree+1) :: Vander, Vander_inv
!!$Integer(int32) :: i, inv_err
!!$
!!$!- End of header -------------------------------------------------------------
!!$
!!$If ( degree .ge. n ) Then
!!$  Write(*,*) 'Error in polyn_fit: degree must be < n',degree,n
!!$Endif
!!$
!!$Vander(:,degree+1) = 1
!!$Do i = degree,1,-1
!!$  Vander(1:n,i) = x*Vander(1:n,i+1)
!!$Enddo
!!$
!!$write(*,'(f8.3)') Vander
!!$
!!$Call Inversion_LU(Vander,Vander_inv,Size(Vander,1),inv_err)
!!$
!!$write(*,'(f8.3)') Vander_inv
!!$
!!$stop
!!$
!!$Vander_inv = Transpose(Vander_inv)
!!$p = Matmul(y,Vander_inv)
!!$End Function polyn_fit

!-----------------------------------------------------------------------------
!+ Calculates the inverse of a square matrix
!-----------------------------------------------------------------------------
Subroutine FINDInv(matrix, inverse, n, errorflag)
!
! Description:
!   !Subroutine to find the inverse of a square matrix
!   !Author : Louisda16th a.k.a Ashwith J. Rego
!   !Reference : Algorithm has been well explained in:
!   !http://math.uww.edu/~mcfarlat/inverse.htm ! ! 
!   !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
!   modified by JL 7/2009 for kind precision
!   "            " 8/31/10 to for free form and match PENTA formatting

Use kind_mod, Only: int32, real64

Implicit None
  
!Declarations
Integer(int32), Intent(In) :: n
Integer(int32), Intent(Out) :: errorflag  !Return error status. -1 for error, 0 for normal
Real(real64), Intent(In), Dimension(n,n) :: matrix  !Input matrix
Real(real64), Intent(Out), Dimension(n,n) :: inverse !Inverted matrix
    
Logical :: FLAG = .TRUE.
Integer(int32) :: i, j, k
Real(real64) :: m
Real(real64), Dimension(n,2*n) :: augmatrix !augmented matrix
!- End of header -------------------------------------------------------------

! Augment input matrix with an identity matrix
Do i = 1, n
  Do j = 1, 2*n
  If (j <= n ) Then
    augmatrix(i,j) = matrix(i,j)
  ElseIf ((i+n) == j) Then
    augmatrix(i,j) = 1._real64
  Else
    augmatrix(i,j) = 0._real64
  Endif
  EndDo
EndDo
    
! Reduce augmented matrix to upper traingular form
Do k =1, n-1
  If (augmatrix(k,k) == 0._real64) Then
    FLAG = .FALSE.
    Do i = k+1, n
      If (augmatrix(i,k) /= 0._real64) Then
        Do j = 1,2*n
          augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
        EndDo
        FLAG = .TRUE.
        Exit
      Endif
      If (FLAG .EQV. .FALSE.) Then
        Write(*,*) "Matrix is non - invertible"
        inverse = 0._real64
        errorflag = -1_int32
        Return
      EndIf
    EndDo
  EndIf
  Do j = k+1, n 
    m = augmatrix(j,k)/augmatrix(k,k)
    Do i = k, 2*n
      augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
    EndDo
  EndDo
EndDo
    
! Test for invertibility
Do i = 1, n
  If (augmatrix(i,i) == 0) Then
    Write(*,*) "Matrix is non - invertible"
    inverse = 0._real64
    errorflag = -1_int32
    Return
  EndIf
EndDo
    
! Make diagonal elements as 1
Do i = 1 , n
  m = augmatrix(i,i)
  Do j = i , (2 * n)! !     
    augmatrix(i,j) = (augmatrix(i,j) / m)
  EndDo
EndDo
    
! Reduced right side half of augmented matrix to identity matrix
Do k = n-1, 1, -1
  Do i =1, k
    m = augmatrix(i,k+1)
    Do j = k, (2*n)
      augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
    EndDo
  EndDo
EndDo
    
! Store answer
Do i =1, n
  Do j = 1, n
    inverse(i,j) = augmatrix(i,j+n)
  EndDo
EndDo
errorflag = 0_int32
End Subroutine FINDinv


!-----------------------------------------------------------------------------
!+ Calculates the inverse of a square matrix
!-----------------------------------------------------------------------------
Subroutine Inversion_LU(A,Y,n,err)
! This subroutine performs the LU inversion of a square matrix.  
! Based on a program by J-P Moreau, Paris:
!   http://jean-pierre.moreau.pagesperso-orange.fr/fortran.html
! Many modifications by JL:
! converted for kind spec, added arguments instead of read, removed allocation
! no implicit variables
! QQ add error flag

Use kind_mod, Only: int32, real64

Implicit None

Integer(int32), Intent(in) :: n

Real(real64) ::  A(n,n)   !real matrix (n x n)
Real(real64) ::  A1(n,n)  !copy of matrix A
Real(real64) ::  Y(n,n)   !real matrix (n x n)
Integer(int32) ::  INDX(n)  !integer vector (n)
Integer(int32) ::  d, rc, i, j, err
!- End of header -------------------------------------------------------------

! Initialize Error flag
err = 0    

! Generate identity matrix
Y = 0.d0
Do i=1, n 
  Y(i,i) = 1.d0 
EndDo

! Save intiial matrix
A1=A

!call LU decomposition routine (only once)
Call LUDCMP(A,n,INDX,D,rc)

!call solver if previous return code is ok
!to obtain inverse of A one column at a time
If ( rc == 0 ) Then
  Do j = 1, n
    Call LUBKSB(A,n,INDX,Y(1,J))
  EndDo
EndIf

!print results or error message
If ( rc == 1 ) Then
  Write(*,*) ' The matrix is singular, no solution !'
EndIf

! Restore initial matrix
A = A1

EnDSubroutine Inversion_LU
           
!-----------------------------------------------------------------------------
!+ Calculates the LU decomposition of a matrix 
!-----------------------------------------------------------------------------
Subroutine LUDCMP(A,N,INDX,D,CODE)
!
! Original code from:
!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
!http://jean-pierre.moreau.pagesperso-orange.fr/fortran.html
!
!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
!
! Modifications by JL to conform to PENTA standards

Use kind_mod, Only: int32, real64

Implicit None

Integer(int32), Parameter :: NMAX = 100_int32
Real(real64), Parameter ::    TINY = 1.5e-18_real64

REAL(real64) ::  AMAX,DUM, SUM, A(N,N),VV(NMAX)
INTEGER(int32) ::  CODE, D, INDX(N), N, imax, k, j, i
!- End of header -------------------------------------------------------------

D = 1 
CODE = 0

Do I=1,N
  AMAX=0.d0
  Do J = 1,N
    If ( DABS(A(I,J)) > AMAX ) Then 
      AMAX=DABS(A(I,J))
    Endif
  EndDo ! j loop
  If ( AMAX < TINY ) Then
    CODE = 1
    Return
  EndIf
  VV(I) = 1._real64 / AMAX
EndDo ! i loop

Do J = 1,N
  Do I = 1,J-1
    SUM = A(I,J)
    Do K=1,I-1
      SUM = SUM - A(I,K)*A(K,J) 
    EndDo ! k loop
    A(I,J) = SUM
  EndDo ! i loop
  AMAX = 0._real64
  Do I = J,N
    SUM = A(I,J)
    Do K=1,J-1
      SUM = SUM - A(I,K)*A(K,J) 
    EndDo ! k loop
    A(I,J) = SUM
    DUM = VV(I)*DABS(SUM)
    If (DUM >= AMAX) Then
      IMAX = I
      AMAX = DUM
    EndIf
  EndDo ! i loop  
   
  If ( J /= IMAX ) Then
    Do K = 1,N
      DUM = A(IMAX,K)
      A(IMAX,K) = A(J,K)
      A(J,K) = DUM
    EndDo ! k loop
    D = -D
    VV(IMAX) = VV(J)
  EndIf

  INDX(J) = IMAX
  If (Dabs(A(J,J))  < TINY) Then
    A(J,J) = TINY
  EndIf

  If ( J /= N ) Then
    DUM = 1._real64 / A(J,J)
    Do I = J+1,N
      A(I,J) = A(I,J)*DUM
    EndDo ! i loop
  EndIf
EndDo ! j loop
EndSubroutine LUDCMP

!-----------------------------------------------------------------------------
!+ Calculates the LU decomposition of a matrix 
!-----------------------------------------------------------------------------
Subroutine LUBKSB(A,N,INDX,B)
!
! Original code from:
!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
!
! Modifications by JL to conform to PENTA standards

Use kind_mod, Only: int32, real64

Implicit None
Real(real64) ::  SUM, A(N,N),B(N)
Integer(int32) ::  INDX(N), ii, i, ll, j, n
!- End of header -------------------------------------------------------------

II = 0

Do I = 1,N
  LL = INDX(I)
  SUM = B(LL)
  B(LL) = B(I)
  If ( II /= 0 ) Then
    Do J=II,I-1
      SUM = SUM - A(I,J)*B(J)
    EndDo ! j loop
  ElseIf ( SUM /= 0._real64 ) Then
    II = I
  EndIf
  B(I) = SUM
EndDo ! i loop

Do I = N,1,-1
  SUM = B(I)
  If ( I < N ) Then
    Do J = I+1,N
      SUM = SUM - A(I,J)*B(J)
    EndDo ! j loop
  EndIf
  B(I) = SUM / A(I,I)
EndDo ! i loop

EndSubroutine LUBKSB

End Module g3df_math_routines_mod
!- End of module header ------------------------------------------------------
