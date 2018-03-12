! Programming name : Possuile2DFlow.f90
! Programming purpose : Compute possuile flow using LBM
! Programming author: Yang Yang :-)

!-----------------------------------------------------------------------------------------------------------
! Simulation Parameter Bolock
MODULE ParametersBlock
implicit none
! Module name: ParametersBlock
! Modul purpose: Specify and compute basic parameters

! Parallel parameter
integer, parameter:: CPUS = 6
  
! Discrete parameters
! Mesh parameters: 
integer, parameter:: Nx = 512, Ny = 512    ! Mesh size 
integer, parameter:: t_max = 1000000       ! Max time step
! Lattice model parameters:
integer, parameter:: D = 2        ! 2 Dimension model
integer, parameter:: Q = 9        ! D2Q9 velocity model
integer, parameter:: Qm = 5       ! D2Q5 velocity model

real*8, parameter:: c = 1.d0      ! lattice speed
real*8, parameter:: delta_x = 1.d0            ! mesh step length
real*8, parameter:: delta_y = delta_x         ! mesh step length
real*8, parameter:: delta_t = delta_x / c     ! time step length

real*8, parameter:: omega_alpha(1:Q) = (/ 4.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0,&
& 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0 /)    ! Weight of distrbution function 
integer, parameter:: e_alpha(1:2,1:Q) = (/0,0, 1,0, 0,1, -1,0, 0,-1, 1,1, -1,1, -1,-1, 1,-1 /)  ! Discrete velocity model

real*8, parameter:: CQ = 1.d0/2.d0                     ! CQ parameter in Mass LBE
real*8, parameter:: J0 = 0.9995d0                      ! J0 parameter in Mass LBE
real*8, parameter:: J_alpha(1:Qm) = (/J0, (1.d0-J0)/4.d0, (1.d0-J0)/4.d0, (1.d0-J0)/4.d0, (1.d0-J0)/4.d0/) ! J_alpha in Mass LBE


! Physics parameters
! Flow parameters
real*8, parameter:: rho_0 = 1.d0    ! Reference Density
real*8, parameter:: Ma = 0.001d0    ! Mach number (Inlet Velocity)
real*8, parameter:: Re = 200        ! Renoyld number (Viscosity)
! Mass Transfer parameters
real*8, parameter:: C0 = 1.d0       ! Inlet concentration
real*8, parameter:: Sc = 2000       ! Schmit number (Diffusion coefficient)
real*8, parameter:: Pe = Re*Sc      ! Peclect number (Diffusion coefficient)

! Computing parameters
! Fluid Flow
real*8, parameter:: Ly = Ny*delta_y                   ! Characteristic length
real*8, parameter:: U0 = 1.d0/sqrt(3.d0)*c*Ma         ! Characteristic Velocity
real*8, parameter:: nu = U0*Ly / Re                   ! Viscosity of fluid 
real*8, parameter:: tau_f = (3.d0*nu) / (delta_t*c**2) + 1.d0/2.d0     ! dimensionless relaxion time of flow
! Mass Transfer
real*8, parameter:: Dm = nu / Sc                   ! Diffusion coefficient
real*8, parameter:: tau_m = Dm / (delta_t*c**2*CQ*(1-J0)) + 1.d0/2.d0  ! dimensionless relaxion time of Mass transfer

! Central parameters:  Nx,Ny,t_max;D,Q,Qm;c,delta_x,delta_y,delta_t,omega_alpha,e_alpha,CQ,J0,J_alpha
! U0,C0,nu,Dm,tau_f,tau_m
END MODULE ParametersBlock
!------------------------------------------------


!------------------------------------------------
! main program
PROGRAM main
use ParametersBlock
implicit none
integer, parameter:: Output_Interval = 10000, Report_Interval = 1000
real*8, parameter:: BreakCriticalCondition = 1.d-5


integer::t=0
real*8 f(Nx+1,Ny+1,Q), f_eq(Nx+1,Ny+1,Q), f_Temp(Nx+1,Ny+1,Q)
real*8 rho(Nx+1,Ny+1), u(Nx+1,Ny+1,D), u_temp(Nx+1,Ny+1,D)
real*8 epsilon

!%% Output Simulation parameters
call OutputParameters()
Write(*,*) 'Programing Paused, Continue?...'
pause

!%% Initial field
call InitialField(f,f_eq,rho,u)
!%% Evolution field main loop
do t = 1,t_max
    call ComputeF(f,f_eq,f_Temp)
    call MacroVars(f,f_Temp,rho,u,u_temp)
    call ComputeF_eq(rho,u,f_eq)
    Call BCs(f,f_eq,rho,u)
    if(mod(t,Output_Interval).EQ.0) then
       Call Outputs(rho,u,t)
    end if
    if(mod(t,Report_Interval).EQ.0 ) then
       Call Erro(u,u_temp,epsilon)
       print 300, t,epsilon
       if(epsilon<BreakCriticalCondition .AND. t>5000) then          
          Print*, 'epsilon < BreakCriticalCondition','Step = ',t
          Exit
       end if
    end if

end do

300 format(1x,'Current Time Step = ',I10,', Current Residuls = ',ES25.15,/)
!%% Report and output computing result 
call SimulationReport(t,epsilon,u)

!%% Output Velocity field and time step as meso file
call MesoField(u,t)


END PROGRAM main
!--------------------------------------------------------------



!--------------------------------------------------------------
!subroutine OutputParameters
SUBROUTINE OutputParameters() 
use ParametersBlock
implicit none
integer i
character(len=20), parameter:: FileName1='SimulationParameters.txt'
!FileName2='MesoParameters.dat'

! Print 
Write(*,100) Nx,Ny,t_max
100 format(1x,'Nx = 'I10,', Ny = ',I10,', t_max = ',I10,/)
Write(*,101) D,Q,Qm
101 format(1x,'D = ',I3,', Q = ',I3,', Qm = ',I3,/)
Write(*,102) c,delta_x,delta_y,delta_t
102 format(1x,'c = ',ES15.5,', delta_x = ',ES15.5,', delta_y = ',ES15.5,', delta_t = ',ES15.5,/)
Write(*,103) omega_alpha
103 format(1x,'Omega_alpha = ',9F10.5,/)
Write(*,*) 'e_alpha = '
do i = 1,D
  Write(*,104) e_alpha(i,:)
  104 format(1x,9I3)
end do

Write(*,105) CQ,J0
105 format(1x,'CQ = ',F10.5,', J0 = 'F10.5,/)
Write(*,106) J_alpha
106 format(1x,'J_alpha ='5F10.5)

Write(*,107) U0,Rho_0,C0
107 format(1x,'U0 = ',F10.5,', Rho_0 = ',F10.5,', C0 = 'F10.5,/)
Write(*,108) nu,Dm
108 format(1x,'nu = ',ES15.6,', Dm = ',ES15.6,/)
Write(*,109) tau_f,tau_m
109 format(1x,'tau_f = ',ES15.6,', tau_m = ',ES15.6,/)

Write(*,110) Re,Sc,Pe
110 format(1x,'Re = ',F12.4,', Sc = ',F12.4,', Pe = ',F12.4,/)

! Output to ReportFile
Open(unit=20,file=Trim(FileName1))
Write(20,100) Nx,Ny,t_max
Write(20,101) D,Q,Qm
Write(20,102) c,delta_x,delta_y,delta_t
Write(20,103) omega_alpha
Write(20,*) 'e_alpha = '
do i = 1,D
  Write(20,104) e_alpha(i,:)
end do

Write(20,105) CQ,J0
Write(20,106) J_alpha

Write(20,107) U0,Rho_0,C0
Write(20,108) nu,Dm
Write(20,109) tau_f,tau_m
Write(20,110) Re,Sc,Pe
Close(20)

! Output Meso parameter file
!Open(unit=30,file=Trim(FileName2))
!Write(30,*) Nx,Ny,t_max,D,Qm,c,delta_x,delta_y,delta_t,e_alpha,CQ,J_alpha,C0,tau_m
!Close(30)

END SUBROUTINE OutputParameters
!---------------------------------------------------------



!---------------------------------------------------------
SUBROUTINE InitialField(f,f_eq,rho,u)
use ParametersBlock
use omp_lib
implicit none
real*8, INTENT(INOUT):: f(Nx+1,Ny+1,Q), f_eq(Nx+1,Ny+1,Q)
real*8, INTENT(INOUT):: rho(Nx+1,Ny+1), u(Nx+1,Ny+1,D)
integer i,j,k
real*8 u_square,e_temp,x

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(u,rho) NUM_THREADS(CPUS)
! Macro Varibles 
do i=1,Nx+1
   do j = 1,Ny+1
     u(i,j,:) = 0
     rho(i,j) = Rho_0
   end do
end do
!$OMP END PARALLEL DO

j = Ny+1
!$OMP PARALLEL DO PRIVATE(i,x) FIRSTPRIVATE(j) SHARED(u) NUM_THREADS(CPUS)
do i = 1,Nx+1
    x = (i-1.d0)/Nx
    u(i,j,1) = U0*16.d0*(x**2)*(1.d0-x)**2
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i,j,u_square,e_temp) SHARED(f,f_eq,u,rho) NUM_THREADS(CPUS)
! Micro Variles
do i=1,Nx+1
   do j = 1,Ny+1
     u_square = u(i,j,1)**2 + u(i,j,2)**2
     do k = 1,Q
        e_temp = u(i,j,1)*dble(e_alpha(1,k)) + u(i,j,2)*dble(e_alpha(2,k)) 
        f_eq(i,j,k) = rho(i,j)*Omega_alpha(k)*( 1.d0 + e_temp*3.d0/c**2 + e_temp**2*9.d0/2.d0/c**4 - u_square*3.d0/2.d0/c**2 )
        f(i,j,k) =  f_eq(i,j,k)
     end do
   end do
end do
!$OMP END PARALLEL DO

END SUBROUTINE InitialField
!--------------------------------------------------------




!--------------------------------------------------------
SUBROUTINE ComputeF(f,f_eq,f_Temp)
use ParametersBlock
use omp_lib
implicit none

real*8, INTENT(INOUT):: f_Temp(Nx+1,Ny+1,Q)
real*8, INTENT(IN):: f(Nx+1,Ny+1,Q), f_eq(Nx+1,Ny+1,Q)

integer i,j,k
integer i_temp,j_temp

!$OMP PARALLEL DO PRIVATE(i,j,k,i_temp,j_temp) SHARED(f,f_eq,F_Temp) NUM_THREADS(CPUS)
do i = 2,Nx
   do j = 2,Ny
     do k = 1,Q
         i_temp = i - e_alpha(1,k)
         j_temp = j - e_alpha(2,k)
         F_Temp(i,j,k) = f(i_temp,j_temp,k) + (f_eq(i_temp,j_temp,k) - f(i_temp,j_temp,k))/tau_f 
     end do
   end do
end do
!$OMP END PARALLEL DO


END SUBROUTINE ComputeF
!--------------------------------------------------------




!--------------------------------------------------------
SUBROUTINE MacroVars(f,f_Temp,rho,u,u_temp)
use ParametersBlock
use omp_lib
implicit none

real*8, INTENT(IN):: f_Temp(Nx+1,Ny+1,Q)
REAL*8, INTENT(OUT):: f(Nx+1,Ny+1,Q), rho(Nx+1,Ny+1),u(Nx+1,Ny+1,D),u_temp(Nx+1,Ny+1,D)
integer i,j,k

!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(u_temp,rho,u,f,f_Temp) NUM_THREADS(CPUS)
! Initial and sum all direction
do i = 2,Nx
   do j = 2,Ny
     u_temp(i,j,:) = u(i,j,:)
     rho(i,j) = 0.d0
     u(i,j,:) = 0.d0
     do k = 1,Q
       f(i,j,k) = f_Temp(i,j,k)
       rho(i,j) = rho(i,j) + f(i,j,k)
       u(i,j,1) = u(i,j,1) + e_alpha(1,k)*f(i,j,k)
       u(i,j,2) = u(i,j,2) + e_alpha(2,k)*f(i,j,k)
     end do
     u(i,j,:) = u(i,j,:) / Rho(i,j)
   end do
end do
!$OMP END PARALLEL DO

END SUBROUTINE 
!--------------------------------------------------------




!--------------------------------------------------------
SUBROUTINE ComputeF_eq(rho,u,f_eq)
use ParametersBlock
use omp_lib
implicit none

real*8, INTENT(IN):: rho(Nx+1,Ny+1),u(Nx+1,Ny+1,D)
real*8, INTENT(OUT):: f_eq(Nx+1,Ny+1,Q)
integer i,j,k
real*8 u_square,e_temp

!$OMP PARALLEL DO PRIVATE(i,j,k,u_square,e_temp) SHARED(u,rho,f_eq) NUM_THREADS(CPUS)
do i = 2,Nx
   do j = 2,Ny
     u_square = u(i,j,1)**2 + u(i,j,2)**2 
     do k = 1,Q
        e_temp = u(i,j,1)*dble(e_alpha(1,k)) + u(i,j,2)*dble(e_alpha(2,k)) 
        f_eq(i,j,k) = rho(i,j)*Omega_alpha(k)*( 1.d0 + e_temp*3.d0/c**2 + e_temp**2*9.d0/2.d0/c**4 - u_square*3.d0/2.d0/c**2 )
     end do
    end do
end do
!$OMP END PARALLEL DO

END SUBROUTINE ComputeF_eq
!--------------------------------------------------------


!--------------------------------------------------------
SUBROUTINE BCs(f,f_eq,rho,u)
use ParametersBlock
use omp_lib
implicit none

real*8, INTENT(INOUT):: f(Nx+1,Ny+1,Q),f_eq(Nx+1,Ny+1,Q),rho(Nx+1,Ny+1),u(Nx+1,Ny+1,D)
integer i,j,k
real*8 u_square,e_temp,x


!%% Inlet non-equilibrium scheme
i = 1
!$OMP PARALLEL DO PRIVATE(j,k,u_square,e_temp) FIRSTPRIVATE(i) SHARED(rho,u,f_eq,f) NUM_THREADS(CPUS)
do j = 2,Ny
   rho(i,j) = rho(i+1,j)
   u(i,j,1) = 0
   u(i,j,2) = 0
   u_square = u(i,j,1)**2 + u(i,j,2)**2 
   do k = 1,Q
     e_temp = u(i,j,1)*dble(e_alpha(1,k)) + u(i,j,2)*dble(e_alpha(2,k))
     f_eq(i,j,k) = rho(i,j)*Omega_alpha(k)*( 1.d0 + e_temp*3.d0/c**2 + e_temp**2*9.d0/2.d0/c**4 - u_square*3.d0/2.d0/c**2 )
     f(i,j,k) = f_eq(i,j,k) + (f(i+1,j,k) - f_eq(i+1,j,k) )
   end do
end do
!$OMP END PARALLEL DO


!%% Outlet fully developed 
i = Nx+1
!$OMP PARALLEL DO PRIVATE(j,k,u_square,e_temp) FIRSTPRIVATE(i) SHARED(rho,u,f_eq,f) NUM_THREADS(CPUS)
do j = 2,Ny
    u(i,j,1) = 0
    u(i,j,2) = 0
    rho(i,j) = rho(i-1,j)
    do k = 1,Q
       e_temp = u(i,j,1)*dble(e_alpha(1,k)) + u(i,j,2)*dble(e_alpha(2,k))
       f_eq(i,j,k) = rho(i,j)*Omega_alpha(k)*( 1.d0 + e_temp*3.d0/c**2 + e_temp**2*9.d0/2.d0/c**4 - u_square*3.d0/2.d0/c**2 )
       f(i,j,k) = f_eq(i,j,k) + (f(i-1,j,k) - f_eq(i-1,j,k) )
    end do 
end do
!$OMP END PARALLEL DO 

!%% Walls non-equilibrium scheme
! %% Bottom 
j = 1
!$OMP PARALLEL DO PRIVATE(i,k,u_square,e_temp) FIRSTPRIVATE(j) SHARED(rho,u,f,f_eq) NUM_THREADS(CPUS)
do i = 1,Nx + 1
   rho(i,j) = rho(i,j+1)
   u(i,j,1) = 0
   u(i,j,2) = 0
   u_square = u(i,j,1)**2 + u(i,j,2)**2 
   do k = 1,Q
     e_temp = u(i,j,1)*dble(e_alpha(1,k)) + u(i,j,2)*dble(e_alpha(2,k))
     f_eq(i,j,k) = rho(i,j)*Omega_alpha(k)*( 1.d0 + e_temp*3.d0/c**2 + e_temp**2*9.d0/2.d0/c**4 - u_square*3.d0/2.d0/c**2 )
     f(i,j,k) = f_eq(i,j,k) + (f(i,j+1,k) - f_eq(i,j+1,k) )
   end do
end do
!$OMP END PARALLEL DO


!%% TOP
j = Ny+1
!$OMP PARALLEL DO PRIVATE(i,k,u_square,e_temp,x) FIRSTPRIVATE(j) SHARED(rho,u,f,f_eq) NUM_THREADS(CPUS)
do i = 1,Nx + 1
   x = (i-1.d0)/Nx
   rho(i,j) = rho(i,j-1)
   u(i,j,1) = U0*16.d0*(x**2)*(1.d0-x)**2
   u(i,j,2) = 0
   u_square = u(i,j,1)**2 + u(i,j,2)**2 
   do k = 1,Q
     e_temp = u(i,j,1)*dble(e_alpha(1,k)) + u(i,j,2)*dble(e_alpha(2,k))
     f_eq(i,j,k) = rho(i,j)*Omega_alpha(k)*( 1.d0 + e_temp*3.d0/c**2 + e_temp**2*9.d0/2.d0/c**4 - u_square*3.d0/2.d0/c**2 )
     f(i,j,k) = f_eq(i,j,k) + (f(i,j-1,k) - f_eq(i,j-1,k) )
   end do
end do
!$OMP END PARALLEL DO

END SUBROUTINE BCs
!--------------------------------------------------------


!--------------------------------------------------------
SUBROUTINE Outputs(rho,u,t)
use ParametersBlock
implicit none

integer, INTENT(IN):: t
real*8, INTENT(IN):: rho(Nx+1,Ny+1),u(Nx+1,Ny+1,D)

real*8 rho_out(Nx+1,Ny+1),u_out(Nx+1,Ny+1,D)
character(len=100) Char_Temp,Nx_Temp,Ny_Temp
integer i,j

Write(Char_Temp,*) t
Write(Nx_Temp,*) Nx+1
Write(Ny_Temp,*) Ny+1
open(20,file='FlowDataAtStep'//Trim(Char_Temp)//'.dat')

Write(20,*) 'VARIABLES = X, Y, U, V, Rho'
Write(20,*) 'Zone T="", I ='//Trim(Nx_Temp)//', J ='//Trim(Ny_Temp)//', F = Point'

! Dimensionless
do i = 1,Nx+1
   do j = 1,Ny+1
      u_out(i,j,:) = u(i,j,:) / U0
      rho_out(i,j) = rho(i,j) / Rho_0
   end do
end do

! Output
do j = 1,Ny+1
   do i = 1,Nx+1
     Write(20,200) dble(i-1)/dble(Ly),dble(j-1)/dble(Ly),u_out(i,j,1),u_out(i,j,2),rho_out(i,j)
   end do
end do
200 format(1x,5ES24.15)

close(20)

END SUBROUTINE Outputs
!--------------------------------------------------------



!--------------------------------------------------------
SUBROUTINE Erro(u,u_temp,epsilon)
use ParametersBlock
use omp_lib
implicit none

real*8, INTENT(IN):: u(Nx+1,Ny+1,D),u_temp(Nx+1,Ny+1,D)
real*8, INTENT(INOUT):: epsilon
integer i,j
real*8:: num = 0.d0, den = 0.d0

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(u,u_temp) REDUCTION(+:num,den) NUM_THREADS(CPUS)
do i = 2,Nx
   do j = 2,Ny
    num = num + (u(i,j,1)-u_temp(i,j,1))**2 + (u(i,j,2)-u_temp(i,j,2))**2
    den = den + u_temp(i,j,1)**2 + u_temp(i,j,2)**2
   end do
end do
!$OMP END PARALLEL DO

epsilon = sqrt(num) / sqrt(den)

END SUBROUTINE Erro
!--------------------------------------------------------



!--------------------------------------------------------
SUBROUTINE SimulationReport(t,epsilon,u)
use ParametersBlock 
implicit none

real*8, INTENT(IN):: t,epsilon,u(Nx+1,Ny+1,D)
integer j,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10

! %% Compute Log File
open(20,file='Report.log')
Write(20,200) t,epsilon
200 format(1x,'Final Time Step:',I10,', Final Residual:',ES25.15)
close(20)

! %% Velocity Profile 
i1 = NINT(0.d0*Nx) +1
i2 = NINT(1.d0/20.d0*Nx) +1
i3 = NINT(1.d0/10.d0*Nx) +1
i4 = NINT(3.d0/20.d0*Nx) +1
i5 = NINT(1.d0/5.d0*Nx) +1
i6 = NINT(3.d0/10.d0*Nx) +1
i7 = NINT(2.d0/5.d0*Nx) +1
i8 = NINT(3.d0/5.d0*Nx) +1
i9 = NINT(4.d0/5.d0*Nx) +1
i10 = Nx+1    
open(30,file = 'VelocityProfile.txt')
Write(30,*) 'Variables = Y,U(X1:X10)'
do j = 1,Ny+1
    Write(30,300) dble(j-1)/dble(Ly) ,u(i1,j,1)/U0,u(i2,j,1)/U0,u(i3,j,1)/U0,u(i4,j,1)/U0,u(i5,j,1)/U0& 
              &,u(i6,j,1)/U0,u(i7,j,1)/U0,u(i8,j,1)/U0,u(i9,j,1)/U0, u(i10,j,1)/U0
end do
300 format(1x,11ES18.10)
close(30)

END SUBROUTINE SimulationReport
!--------------------------------------------------------


!--------------------------------------------------------
!: This Subroutine is really important, play a role of bridge
SUBROUTINE Mesofield(u,t)
use ParametersBlock
implicit none

real*8, INTENT(IN):: u(Nx+1,Ny+1,D)
integer, INTENT(IN):: t
open(100,file = 'MesoField.dat')
Write(100,*) u,t

close(100)
END SUBROUTINE Mesofield
!--------------------------------------------------------
