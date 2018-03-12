module CONJ_GRAD
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-5
!-----------------------------------------------------
!  Description : 共轭梯度法
!    
!-----------------------------------------------------
!  Parameters  :
!      1.    IMAX--最大允许迭代次数  
!      2.    tol--误差容限
!  
!  Contains    :
!      1.    solve 迭代法方法函数
!      2.    
!      3.    dr(r,N) 计算向量长度平方函数
!      4.     Ar(A,r,N)  计算矩阵乘以向量函数，返回向量
!      5.    rAr(A,r,N)   计算（Ar,r）函数
!-----------------------------------------------------
!  Post Script :
!      1.    里面的三个函数，可以简化程序，同时可以用在其他地方
!      2. 
!-----------------------------------------------------

implicit real*8(a-z)
integer::IMAX = 5000
real*8::tol = 1d-16


contains

subroutine solve(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  共轭梯度法
!               用于计算方程 AX=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.  A,b 意义即  AX=b
!       2.  x0迭代初值
!       3.  N 方程的维数
!  Output parameters  :
!       1. x 方程的解
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::N
integer::i,j,k

real*8::A(N,N),b(N),x(N),x0(N)
real*8::r0(N),r1(N),p0(N),p1(N)

real*8::x1(N),x2(N)

! Message
print*, 'Enter Subroutine solve...'
!-------------------------------------------------------------------------------------------------------

! Open grid file
open(102,FILE='Solver.dat')
!-------------------------------------------------------------------------------------------------------
!写入标题
  write(102,501)
  501 format(1x,//,18x,'共轭梯度法中间结果',//)
!-------------------------------------------------------------------------------------------------------

! Initializing value for iterative
do i = 1,N
   x0(i) = 0.d0
end do

!--------------------------------------------------------------------------------------------------------

x1=x0

r0=b-Ar(A,x1,N)

p0=r0

do k=1,IMAX
    
    
    tmp1=dr(r0,N)
    tmp2=rAr(A,p0,N)
    
    afa=tmp1/tmp2
    
    x2=x1+afa*p0
    
    

    
    
    !如果r0接近于0，则停止迭代
    !该部分算法在《数值分析原理》（李庆扬、关治、
    !白峰彬编）中叙述较为详细
    dr_s=dsqrt(dr(r0,N))
    
    if (dr_s<tol) exit
    
    
    r1=r0-afa*Ar(A,p0,N)
    
    tmp3=dr(r1,N)
    
    beta=tmp3/tmp1
    
    p1=r1+beta*p0
    
    ! Residual in interative
    write(102,502) k, tmp3
    502 format(1x,'k = ',1I8,',',5x,'residual = ',ES15.8)
    
    !全部更新
    r0=r1
    p0=p1
    x1=x2
         
end do

     x=x2
 
! Message
print*, 'Subroutine solve exit...'
print*, ''
!-------------------------------------------------------------------------------------------------------

! Subroutine end    
end subroutine solve




function dr(r,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  计算向量长度平方  (r,r)  
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     r向量
!       2.     N维数
!  Output parameters  :
!       1.     dr 长度平方
!       2.
!  Common parameters  :
!
!----------------------------------------------------
                                                                  
implicit real*8(a-z)
integer::N,i
real*8::r(N),dr

s=0
do i=1,N
  s=s+r(i)**2
end do
dr=s
end function dr

function Ar(A,r,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  !计算  A*r,返回 N维向量    
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     r向量
!       2.     N维数
!       3.     A矩阵
!  Output parameters  :
!       1.     Ar返回向量
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)
integer::i,N
real*8::A(N,N),r(N),temp(N),Ar(N)

temp=0

do i=1,N
  do j=1,n
    temp(i)=temp(i)+A(i,j)*r(j)
  end do
end do
Ar=temp
end function ar

function v1v2(v1,v2,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-7-29
!-----------------------------------------------------
!  Purpose     : 向量点乘   v1v2=v1(1)*v2(1)+v1(2)*v(2)+...
! 
!  Post Script :
!       1.
!       2.
!       3.
!
!-----------------------------------------------------
!  Input  parameters  :
!       1.   v1,v2 向量
!       2.   N 向量维数
!  Output parameters  :
!       1.   v1,v2向量点乘值
!       2.
!  Common parameters  :
!       1.
!       2.
!----------------------------------------------------

implicit real*8(a-z)
integer::n

real*8::v1(n),v2(n)

integer::i

v1v2=0
do i=1,n

 v1v2=v1v2+v1(i)*v2(i)

end do
end function




function rAr(A,r,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  !计算（Ar,r）,返回标量  
!    
!-----------------------------------------------------
!  Input  parameters  :
!       1.     r向量
!       2.     N维数
!       3.     A矩阵
!  Output parameters  :
!       1.     Ar返回标量
!       2.
!  Common parameters  :
!----------------------------------------------------

implicit real*8(a-z)
integer::i,N
real*8::A(N,N),r(N),temp(N)

temp=Ar(A,r,N)
rAr=v1v2(r,temp,N)
end function rAr


end module CONJ_GRAD
