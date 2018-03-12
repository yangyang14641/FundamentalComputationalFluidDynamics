MODULE READ_GRID
! Read Grid data
CONTAINS

SUBROUTINE ReadGrid(NP,NE,X,Y,NOD)
! This Subroutine is to read mesh in Finite Element Method
implicit none 
! variables
integer NP,NE 
real*8, allocatable,INTENT(INOUT):: X(:),Y(:)
integer, allocatable,INTENT(INOUT)::  NOD(:,:)

! local Variables
integer i
integer status_msg
!-------------------------------------------------------------------------------------------------------


! Message
print*, 'Enter Subroutine: ReadGrid...'
!-------------------------------------------------------------------------------------------------------

! Open grid file
open(10,FILE='GRID.DAT')
!-------------------------------------------------------------------------------------------------------

! read nodes number
read(10,1001) NP,NE
print*, 'NUmber of Points =',NP,', Number of Element = ',NE
!-------------------------------------------------------------------------------------------------------

! Alloc memory 
ALLOCATE(X(1:NP), STAT = status_msg )
if (status_msg == 0) then
  print*, 'Allocate X(NP) success!'
else
  print*, 'Allocate X(NP) fail!'
end if
ALLOCATE(Y(1:NP), STAT = status_msg)
if (status_msg == 0) then
  print*, 'Allocate Y(NP) success!'
else
  print*, 'Allocate Y(NP) fail!'
end if
ALLOCATE(NOD(1:3,1:NE), STAT = status_msg)
if (status_msg == 0) then
  print*, 'Allocate NOD(1:3,NE) success!'
else
  print*, 'Allocate NOD(1:3,NE) fail!'
end if
!-------------------------------------------------------------------------------------------------------

! Read nodes' position
do i=1,NP
   read(10,1002) X(i),Y(i)
end do

!print*, X,Y
!-------------------------------------------------------------------------------------------------------

! Read elements' topology
do i=1,NE
  read(10,1003) NOD(1,i),NOD(2,i),NOD(3,i)
end do
!-------------------------------------------------------------------------------------------------------

close(10)
1001 format(1X,2I10)
1002 format(1X,E13.6,2X,E13.6)
1003 format(1X,3I10)
!-------------------------------------------------------------------------------------------------------


! Message
print*, 'Subroutine: ReadGrid exit...'
print*, ''
!-------------------------------------------------------------------------------------------------------


! Program end
END SUBROUTINE ReadGrid
END MODULE READ_GRID
