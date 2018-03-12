!SUBROUTINE CHASE.f90
!**************************************************************
!**************************************************************
! Subroutine Name: Chase
! Purpose : Use Thomas Algorithm to solve spares Linear Algebra Equations
! Author: Yang Yang
! October. 2013. 
!**************************************************************
!**************************************************************
      
SUBROUTINE Chase(c,e,b,f,Nx,x)
IMPLICIT NONE
        
INTEGER Nx
REAL*8 e(2:Nx),b(Nx),c(1:Nx-1),f(Nx),x(Nx)
! LU diag
REAL*8 l(2:Nx),u(Nx),y(Nx)
! Control Varibles:
INTEGER i
        
! Thomas Algorithm:
        
  ! step one: LU Decompose
  ! x-->  x
  !       |
  !       V
  !       x --> x
  u(1)= b(1)
  DO i=2,Nx
       l(i) = e(i) / u(i-1)
       u(i) = b(i) - l(i)*c(i-1)
  END DO
             
  ! step two: Backward (pop)
  !   |1           |    |f(1)|
  !   |x 1         |    |f(2)|
  !   |  x 1       |    |f(3)|
  ! L=|    . .     |  f=| .  |
  !   |      . .   |    | .  |
  !   |         x 1|    |f(N)|
  ! so y(i)=[y(1),y(2)......y(N)]
  y(1) = f(1)
  DO i=2,Nx
      y(i) = f(i) - l(i)*y(i-1)
  END DO
             
  !   |x x              |    |y(1)|
  !   |  x x            |    |y(2)|
  !   |    x x          |    |y(3)|
  !   |      . .        |    |  . |
  ! U=|        . .      | y= |  . |
  !   |          . .    |    |  . |
  !   |            . .  |    |  . |
  !   |              x x|    |y(N)|
  x(Nx) = y(Nx) / u(Nx)
  DO i=Nx-1,1,-1
       x(i) = ( y(i)-c(i)*x(i+1) ) / u(i)
  END DO
             
END SUBROUTINE CHASE