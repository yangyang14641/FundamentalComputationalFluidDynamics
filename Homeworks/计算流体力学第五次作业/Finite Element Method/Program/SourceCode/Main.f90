PROGRAM Main
! This is main program to compute potential flow by using Finite Element Method
USE READ_GRID          ! Module read grid 
USE FIND_BOUND         ! Module find boundary
USE CONJ_GRAD          ! Module Conjugate Gradient Solver for symmetry positive define system
implicit none
!------------------------------------------------------------------------

! Mesh variables
integer NP,NE                              ! Number of Nodes and Elements
real*8, allocatable:: X(:),Y(:)            ! Position of Nodes
integer, allocatable::  NOD(:,:)           ! Topology of Elements

! Record boundary nodes
integer:: N_ab = 0
integer, allocatable::  ab(:)
integer:: N_bc = 0
integer, allocatable::  bc(:)
integer:: N_de = 0
integer, allocatable::  de(:)
integer:: N_ea = 0
integer, allocatable::  ea(:)


! local index variables
integer i,j,ii
integer status_msg

! Stiffness Matrix and right column vector
real*8, allocatable:: A(:,:), f(:)

real*8 A_e                    ! Area of element
real*8 b1,b2,b3,c1,c2,c3      ! interpolation coefficient
integer i_e1,i_e2,i_e3        ! index of elements' nodes 

! Punishment factor for Dirichelet boundary condition
real*8,parameter:: punishment = 1D8


! Result variables
real*8, allocatable:: phi(:),u(:),p(:)


!------------------------------------------------------------------------

! Message
print*, 'Program Main exit...'
!------------------------------------------------------------------------

! Read Grid data
call ReadGrid(NP,NE,X,Y,NOD)              ! Call sbroutine to read grid
!------------------------------------------------------------------------

! Find Boundary Nodes 
call FindBoundary(NP,X,Y,N_ab,N_bc,N_de,N_ea,ab,bc,de,ea)
!------------------------------------------------------------------------

! Generate Element Stifness Matrix and Assemble Total Stiffness Matrix
Allocate(A(1:NP,1:NP), STAT = status_msg)
print*, 'memory allocate state for A(:,:) is:("0" means success)  ',status_msg

do i = 1,NP
   do j = 1,NP
      A(i,j) = 0.d0         ! initializing stiffness matrix
   end do
end do

do i = 1,NE
   ! load nodes index
   i_e1 = NOD(1,i)
   i_e2 = NOD(2,i)
   i_e3 = NOD(3,i)
   
   ! computing element's area
   A_e = 5.d-1*( (X(i_e2)-X(i_e1))*(Y(i_e3)-Y(i_e1)) - (Y(i_e2)-Y(i_e1))*(X(i_e3)-X(i_e1)) )
   
   ! computing interpolation coefficient b_{i}
   b1 = 5.d-1 * 1.d0/A_e*( Y(i_e2)-Y(i_e3) )
   b2 = 5.d-1 * 1.d0/A_e*( Y(i_e3)-Y(i_e1) )
   b3 = 5.d-1 * 1.d0/A_e*( Y(i_e1)-Y(i_e2) )

   ! computing interpolation coefficient c_{i}
   c1 = 5.d-1 * 1.d0/A_e*( X(i_e3)-X(i_e2) )
   c2 = 5.d-1 * 1.d0/A_e*( X(i_e1)-X(i_e3) )
   c3 = 5.d-1 * 1.d0/A_e*( X(i_e2)-X(i_e1) )
   
   ! assemble total stiffness matrix
   A(i_e1,i_e1) = A(i_e1,i_e1) + b1**2 + c1**2 
   A(i_e1,i_e2) = A(i_e1,i_e2) + b2*b1 + c2*c1
   A(i_e1,i_e3) = A(i_e1,i_e3) + b3*b1 + c3*c1
     
   A(i_e2,i_e1) = A(i_e2,i_e1) + b1*b2 + c1*c2
   A(i_e2,i_e2) = A(i_e2,i_e2) + b2**2 + c2**2
   A(i_e2,i_e3) = A(i_e2,i_e3) + b3*b2 + c3*c2
     
   A(i_e3,i_e1) = A(i_e3,i_e1) + b1*b3 + c1*c3
   A(i_e3,i_e2) = A(i_e3,i_e2) + b2*b3 + c2*c3
   A(i_e3,i_e3) = A(i_e3,i_e3) + b3**2 + c3**2

end do

!------------------------------------------------------------------------

! Adding Dirichelet boundary condition in Total Stiffness Matrix and creat right column vector
Allocate(f(1:NP), STAT = status_msg)
print*, 'memory allocate state for f(:) is:("0" means success)  ',status_msg

do i = 1,NP       ! Initializing f(:)
    f(i) = 0.d0
end do


do i = 1,N_ab     ! boundary ab
     ii = ab(i)   ! get node's ID in boundary
     A(ii,ii) = A(ii,ii) + punishment
     f(ii) = punishment * 0.d0
end do

do i = 1,N_bc     ! boundary bc
     ii = bc(i)   ! get node's ID in boundary
     A(ii,ii) = A(ii,ii) + punishment
     f(ii) = punishment * 0.d0
end do

do i = 1,N_de     ! boundary de
     ii = de(i)   ! get node's ID in boundary
     A(ii,ii) = A(ii,ii) + punishment
     f(ii) = punishment * 2.d0
end do

do i = 1,N_ea     ! boundary ea
     ii = ea(i)   ! get node's ID in boundary
     A(ii,ii) = A(ii,ii) + punishment
     f(ii) = punishment * Y(ii)
end do

!------------------------------------------------------------------------

! Solve Total Stiffness Matrix
Allocate(phi(1:NP), STAT = status_msg)
print*, 'memory allocate state for phi(:) is:("0" means success)  ',status_msg
Allocate(u(1:NP), STAT = status_msg)
print*, 'memory allocate state for u(:) is:("0" means success)  ',status_msg
Allocate(p(1:NP), STAT = status_msg)
print*, 'memory allocate state for p(:) is:("0" means success)  ',status_msg


call solve(A,f,phi,NP)
!------------------------------------------------------------------------

! Regroup data


!------------------------------------------------------------------------

! Output Result 
call Output(NP,NE,X,Y,Phi,NOD)
!------------------------------------------------------------------------

! Free Memory
DEALLOCATE(X,Y,NOD)
DEALLOCATE(ab,bc,de,ea)
!------------------------------------------------------------------

! Message
print*, 'Program Main exit...'
print*, ''
!------------------------------------------------------------------------

! Program end
END PROGRAM Main
