MODULE FIND_BOUND
! Find Mesh Boundarys
CONTAINS

SUBROUTINE FindBoundary(NP,X,Y,N_ab,N_bc,N_de,N_ea,ab,bc,de,ea)
! This subroutine is to find mesh boundary 
implicit none

! Mesh Information
integer NP
real*8 X(1:NP),Y(1:NP)

! Record boundary nodes
integer N_ab 
integer, allocatable::  ab(:)
integer N_bc 
integer, allocatable::  bc(:)
integer N_de 
integer, allocatable::  de(:)
integer N_ea 
integer, allocatable::  ea(:)

! Local variables
integer:: i=0, i_ab=1, i_bc=1, i_de=1, i_ea=1
integer status_msg
!-------------------------------------------------------------------------------------------------------


! Message
print*, 'Enter Subroutine FindBoundary...'
!-------------------------------------------------------------------------------------------------------

! Find Boundary Nodes
do i = 1,NP
  ! ab
  if( ABS(Y(i) - 0.D0) < 1D-10 ) then
   N_ab = N_ab + 1
  end if

  ! bc
  if( ABS(X(i)**2 + Y(i)**2 - 1) < 1D-4 .and. ABS(Y(i) - 0.D0) > 1D-10) then
   N_bc = N_bc + 1
  end if

  ! de
  if( ABS(Y(i) - 2.D0) < 1D-10 ) then
   N_de = N_de + 1
  end if

  ! ea
  if( ABS(X(i) + 3.5D0) < 1D-10 .and. ABS(Y(i) - 0.D0) > 1D-10 .and. ABS(Y(i) - 2.D0) > 1D-10 ) then
   N_ea = N_ea + 1
  end if

end do

print*, 'N_ab=', N_ab
print*, 'N_bc=', N_bc
print*, 'N_de=', N_de
print*, 'N_ea=', N_ea


ALLOCATE(ab(1:N_ab), STAT = status_msg )
print*, 'memory allocate state for ab is:("0" means success)  ',status_msg
ALLOCATE(bc(1:N_bc), STAT = status_msg )
print*, 'memory allocate state for bc is:("0" means success)  ',status_msg
ALLOCATE(de(1:N_de), STAT = status_msg )
print*, 'memory allocate state for de is:("0" means success)  ',status_msg
ALLOCATE(ea(1:N_ea), STAT = status_msg )
print*, 'memory allocate state for ea is:("0" means success)  ',status_msg

do i = 1,NP
  ! ab
  if( ABS(Y(i) - 0.D0) < 1D-10 ) then
   ab(i_ab) = i
   i_ab = i_ab + 1
  end if

  ! bc
  if( ABS(X(i)**2 + Y(i)**2 - 1) < 1D-4 .and. ABS(Y(i) - 0.D0) > 1D-10) then
   bc(i_bc) = i
   i_bc = i_bc + 1
  end if

  ! de
  if( ABS(Y(i) - 2.D0) < 1D-10 ) then
   de(i_de) = i
   i_de = i_de + 1
  end if

  ! ea
  if( ABS(X(i) + 3.5D0) < 1D-10 .and. ABS(Y(i) - 0.D0) > 1D-10 .and. ABS(Y(i) - 2.D0) > 1D-10 ) then
   ea(i_ea) = i
   i_ea = i_ea + 1
  end if

end do
!------------------------------------------------------------------------------------------------------------------------

! Message
print*, 'Subroutine FindBoundary exit...'
print*, ''
!------------------------------------------------------------------------------------------------------------------------

! Program end
END SUBROUTINE FindBoundary
END MODULE FIND_BOUND
