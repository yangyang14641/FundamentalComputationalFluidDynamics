!  OndedimensionUnsteadyHeatConduct.f90
!
!  Program:
!  OndedimensionUnsteadyHeatConduct - Compute Onedimension Heat Transfer Problem
!  Author: Yang Yang 
!  October. 2013. 

!****************************************************************************
!
!  PROGRAM: OndedimensionUnsteadyHeatConduct
!
!  PURPOSE: Compute Onedimension Heat Conduct Problem in the Wall with Convection Boundary Condition
!
!****************************************************************************

PROGRAM OndedimensionUnsteadyHeatConduct
IMPLICIT NONE

! Basic Physics Parameters: data for metal iron 
REAL*8:: k=80.2,rho=7870,cp=447,L=10,h=300,Tf=300,Ti=900   
! Basic Mesh Parameters:
INTEGER Nx,Nt
PARAMETER(Nx=100,Nt=100)
REAL*8:: dx,dt,x(Nx),t(Nt)
! Temperature Field Variables:
REAL*8,DIMENSION(Nt,Nx)::Temperature 
! Process Varibles:
REAL*8:: e(2:Nx),b(1:Nx),c(1:Nx-1),f(1:Nx),Temp(Nx),Fo,Bi
! Control Varibles:
INTEGER:: i,j

 
! out put file directory
OPEN(2,file='space_time.txt')
OPEN(3,file='result.txt')

! Mesh step length
dx=L/(Nx-1)
dt=2000
x(1) = 0
DO j=2,Nx
   x(j) = x(j-1) + dx
END DO
t(1) = 0
DO j=2,Nt
   t(j) = t(j-1) + dt
END DO

! Compute Fo and Bi (dimensionless)
Fo = (k/(rho*cp)) *dt / (dx**2)
Bi = h*dx / (k)

! Initial Temperature Field
DO i=1,Nt
   DO j=1,Nx
      Temperature(i,j) = 0
   END DO
END DO

DO j=1,Nx
   Temperature(1,j) = Ti
END DO


! main time loop
! |b1 c1                 |
! |e2 b2 c2              |
! |   e3 b3 c3           |
! |      e4 b4 c4        |
! |          .  .  .     |
! |             .  . cn-1|
! |                en  bn|
DO i=2,Nt
   ! Assemble Sparse matrix A
    DO j = 2,Nx-1
       c(j) = -Fo
       e(j) = -Fo
    END DO
    c(1)   = -2*Fo
    e(Nx)  = -2*Fo
    
    DO j=2,Nx-1
       b(j) = (1+2*Fo)
    END DO
    b(1)  = 1+2*Bi*Fo+2*Fo 
    b(Nx) = 1+2*Bi*Fo+2*Fo
    
    DO j=2,Nx-1
       f(j)= Temperature(i-1,j)
    END DO
    f(1)  = 2*Bi*Fo*Tf + Temperature(i-1,1) 
    f(Nx) = 2*Bi*Fo*Tf + Temperature(i-1,Nx)
    DO
    j=1,Nx
      Temp(j) = Temperature(i-1,j)
    END DO
    ! Inner loop to update temperature field:
    CALL chase(c,e,b,f,Nx,Temp)
    ! Save Result
    DO j=1,Nx
      Temperature(i,j) = Temp(j)
    END DO
    ! Print current time step
    PRINT*, "Current Time Step:",i
END DO
 
! save result to file
WRITE(2,200)t(:),x(:)
DO i=1,Nt
   WRITE(3,100)Temperature(i,:)
END DO
100 FORMAT(100(f15.8)) 
200 FORMAT(100(f12.5))

! Body of OndedimensionUnsteadyHeatConduct
PRINT *, 'Computing Job Compelete!'
PAUSE
END PROGRAM OndedimensionUnsteadyHeatConduct

