!Refer the file NM-Exe-4A for the data used
PROGRAM Interpolation2
IMPLICIT NONE

!Declaring Arrays and Variables
REAL,ALLOCATABLE::X(:,:),A(:),F(:),XD(:)
INTEGER::N,i,j,k
REAL::P,D,M,FU
PRINT *,"ENTER THE NUMBER OF POINTS"
READ *,N
read *,P
ALLOCATE(X(N,N),A(N),F(N),XD(N))


!Reading the File
OPEN (23, FILE="data.txt", STATUS="OLD")
DO i=1,N
READ (23,*)XD(i),F(i)
END DO
CLOSE(23)

!Making a lower triangular matrix
X=0
do i=1,n
 do j=1,i
 x(i,j)=1
 end do
end do 


!Forming the division matrix via Newton division method
do i=2,n
 do j=2,i
 x(i,j)=x(i,j-1)*(xd(i)-xd(j-1))
 end do
end do 

!The back substitution
a(1)=f(1)
do i=2,n
d=0
 do j=1,i
 d=d+a(j)*x(i,j)
 end do 
a(i)=1/x(i,i)*(f(i)-d)
end do

!Generation of the final polynomial
fu=a(1)
do i=2,n
m=1
 do j=1,i-1
 m=m*(P-xd(j))
 end do
fu=fu+m*a(i)
end do

print *,fu

END PROGRAM
