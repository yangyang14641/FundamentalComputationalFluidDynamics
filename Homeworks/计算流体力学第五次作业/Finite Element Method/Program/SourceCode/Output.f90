SUBROUTINE Output(NP,NE,X,Y,Phi,NOD)
! This program is to output data
implicit none

integer NP,NE 
real*8  X(1:NP),Y(1:NP),Phi(1:NP)
integer NOD(1:3,1:NE)

! local variables
integer i
!-------------------------------------------------------------------------------------------------------

!Message
print*, 'Enter Subroutine Output...'
!-------------------------------------------------------------------------------------------------------

!Output data in tecplot format
open(10,file='Result.dat')

Write(10,1001) 'TITLE = "FEM:Potential Flow"'
Write(10,1002) 'VARIABLES = "X", "Y", "Phi"'
Write(10,1003) 'ZONE T = "TRIANGLES", NODES =',NP,', ELEMENTS =',NE,', DATAPACKING =POINT, ZONETYPE = FETRIANGLE'
!-------------------------------------------------------------------------------------------------------

! Write nodes' position
do i=1,NP
   write(10,1004) X(i),Y(i),Phi(i)
end do

write(10,*) ''
!-------------------------------------------------------------------------------------------------------

! Write elements' topology
do i=1,NE
  write(10,1005) NOD(1,i),NOD(2,i),NOD(3,i)
end do
!-------------------------------------------------------------------------------------------------------

1001 format(1x,A)
1002 format(1x,A)
1003 format(1x,A,1I6,A,1I6,A)
1004 format(1X,3ES15.6)
1005 format(1X,3I10)

close(10)
!-------------------------------------------------------------------------------------------------------

! Message
print*, 'Subroutine Output exit...'
print*, ''
!-------------------------------------------------------------------------------------------------------

! Subroutine end
END SUBROUTINE OutPut
