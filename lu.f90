PROGRAM ludecomp
IMPLICIT NONE

!Declaring variables and allocating
REAL,ALLOCATABLE::A(:,:),B(:),L(:,:),U(:,:),X(:),Y(:)
REAL::AB
INTEGER::I,J,K,N
PRINT *,"Enter the size of matrix used"
READ *,N
ALLOCATE(A(N,N),B(N),L(N,N),U(N,N),X(N),Y(N))

!Reading the file for the matrix
OPEN(10, FILE="lumat.dat", STATUS="OLD")
DO I=1,N
READ (10,*)A(:,I),B(I)
END DO
CLOSE(10)

!Performing LU Decomposition

!Crout's algorithm
DO I=1,N
L(I,I)=1
END DO

DO I=1,N
 DO J=1,I
 AB=0
  DO K=1,J-1
  AB=AB+L(K,J)*U(I,K)
  END DO
  U(I,J)=A(I,J)-AB
 END DO
 
  DO J=I+1,N
   AB=0
   DO K=1,J-1
   AB=AB+L(K,J)*U(I,K)
   END DO
   L(I,J)=(A(I,J)-AB)/U(I,I)
  END DO
END DO 

!Solving equations
DO I=1,N
 AB=0
 DO J=1,I-1
  AB=AB+L(J,I)*Y(J)
 END DO
 Y(I)=(B(I)-AB)/L(I,I)
END DO

DO I=N,1,-1
 AB=0
 DO J=N,I+1,-1
  AB=AB+U(J,I)*X(J)
 END DO
 X(I)=(Y(I)-AB)/U(I,I)
END DO


!Testing outputs
PRINT *,"================================================================="
PRINT *,"The given matrix is"
DO I=1,N
PRINT *,(A(J,I),J=1,N)
END DO
PRINT *,"The and the constant part is"
PRINT *,B
PRINT *,"================================================================="
PRINT *,"Lower triangular matrix obtained is"
DO I=1,N
PRINT *,(L(J,I),J=1,N)
END DO
PRINT *,"================================================================="
PRINT *,"Upper triangular matrix obtained is"
DO I=1,N
PRINT *,(U(J,I),J=1,N)
END DO
PRINT *,"================================================================="
PRINT *,"The solutions of the equation is"
PRINT *,X
PRINT *,"================================================================="
DEALLOCATE(A,B,L,U,X,Y)
END PROGRAM

! Enter the size of matrix used
!3
! =================================================================
! The given matrix is
!   1.00000000       5.00000000       2.00000000    
!   6.00000000       3.00000000       7.00000000    
!   7.00000000       2.00000000       10.0000000    
! The and the constant part is
!   2.00000000       4.00000000       5.00000000    
! =================================================================
! Lower triangular matrix obtained is
!   1.00000000       0.00000000       0.00000000    
!   6.00000000       1.00000000       0.00000000    
!   7.00000000       1.22222221       1.00000000    
! =================================================================
! Upper triangular matrix obtained is
!   1.00000000       5.00000000       2.00000000    
!   0.00000000      -27.0000000      -5.00000000    
!   0.00000000       0.00000000       2.11111116    
! =================================================================
! The solutions of the equation is
! 0.122807026      0.228070185      0.368420988    
! =================================================================




