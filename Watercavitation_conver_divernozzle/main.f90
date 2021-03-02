!-----------------------------------------------------------------` 
! Reacting Flows
! Project - 3, Water cavitation
! Rohit K S S Vuppala
! rvuppal@ncsu.edu
! 200211761
!-----------------------------------------------------------------

program main
implicit none

real(8) :: u0,t0,P0stag,pratio,Pback,vp0,yv0
real(8),allocatable :: X(:),S(:),Sm(:,:),rho(:),v(:,:),F(:,:),c(:),dudv(:,:,:),&
                      rhs(:,:)
real(8) :: h,y,res(3),res0(3)
integer :: i,j,nx,iflux,niter,iter,flag
real(8) :: Rvap,dt,leng,tol,pin
!----------
!PRIMITIVES ARE:
! Y-vapor, Pressure,Velocity
!CONSERVED :
!rho*Yv,rho,rho*u
!------------------
 
real(8),parameter :: r0 = 1000.0d0  
real(8),parameter :: p0 = 101325.0d0
real(8),parameter :: k0 = 2.15d9
real(8),parameter :: n0=  7.15d0
               
flag = 0
tol  = 1d-6
niter = 10000000
!-------vapour realtion-----------
Rvap = 461.9d0 
vp0  = 3169.0d0
!---------------------------------

nx   = 150
leng = 1.5d0
h    = leng/nx
iflux = 3     !1-Loui Stefan,2-LDFSS,3-Loui Stefan modified,4-LDFSS modified
dt = 1d-6

!---Input values---------------------------------------------

P0stag = 4.0d7
u0     = 10.0d0
t0     = 300.0d0
Pratio = 0.9d0
Pback  = Pratio*P0stag
yv0    = 1d-7
!----------------------------------------------------

allocate(S(0:nx+1))
allocate(X(0:nx+1))
allocate(rho(0:nx+1))
allocate(v(3,0:nx+1))
allocate(Sm(3,1:nx))
allocate(F(3,0:nx) )
allocate(c(0:nx+1) )
allocate(dudv(3,3,1:nx))
allocate(rhs(3,1:nx))

 do i=0,nx+1       !faces are from 0 to nx
    y = i*h
  X(i)= -h/2 + i*h

         S(i) = ((1.0d0+4.0d0*(y-0.5d0)**2))
 
 end do

open(file='grid.dat',unit=1)
do i=1,nx
   
   write(1,*) X(i),S(i)

enddo

close(1)
!--------initialise the values--------------------
call initial
!print*, v(:,1)
!stop 
open(file='convergence.dat',unit=3)
print*,'------------------<<<<<<<<<Initialisation complete>>>>>>>>-----------------'
do iter=1,niter

	call therm
        !print*, c(2)
        call boundary
        call therm
! do i=0,nx+1
!    print*, c(i),v(:,i),rho(i)
! enddo
! stop
!        print*, c(0)
        call flux
        call source
        call update
        call therm
        call residual 
       if(flag.eq.1) exit
       if(mod(iter,1000).eq.0) call solution
enddo
       call solution
contains
!-------------------------------------------------------------
subroutine solution

  open(file='solution.dat',unit=2,status = 'unknown')
  write(2,*) 'variables= "x","Vapour fraction","Pressure","Velocity","Density","Mach","Liquid Fraction"'
  do i=1,nx
     write(2,*) x(i),v(:,i),rho(i),v(3,i)/c(i),1.0-v(1,i)  
  enddo
  close(2)
endsubroutine

!-----------------------------------------------------------
subroutine residual
integer :: i,j
    
    if(iter.eq.1) res0=res 
    if(mod(iter,1000).eq.0) write(*,*) iter, res,abs(sum(res)/sum(res0))
    write(3,*) iter, res,abs(sum(res)/sum(res0))
    if(abs(sum(res)/sum(res0)).le.tol) flag =1
    !stop 
   !call sleep(1)
endsubroutine
!-------------------------------------------------------------

function solver(a,b)
integer,parameter::n = 3
real(8)::a(n,n),s(n),b(n),x(n),r,rmax,smax,sum,z,solver(3)
integer::p(n), pk,i,j,k
!print*, iter
!print*, 'hi'


      do  i=1,n
         p(i) = i
         smax = 0.0d0
         do  j=1,n
            smax = max(smax,dabs(a(i,j)))
         enddo
         s(i) = smax
     enddo

      do  k=1,n-1
         rmax = 0.0d0
         do i=k,n
            r = dabs(a(p(i),k))/s(p(i))
            if (r .gt. rmax) then
               j = i
               rmax = r
            endif
      enddo

         pk = p(j)
         p(j) = p(k)
         p(k) = pk

         do i=k+1,n
            z = a(p(i),k)/a(p(k),k)
            a(p(i),k) = z
            do  j=k+1,n
               a(p(i),j) = a(p(i),j) - z*a(p(k),j)
            enddo
        enddo
      enddo

      do  k=1,n-1
         do  i=k+1,n
            b(p(i)) = b(p(i)) - a(p(i),k)*b(p(k))
         enddo
      enddo
      do i=n,1,-1
         sum = b(p(i))
         do j=i+1,n
            sum = sum - a(p(i),j)*x(j)
    enddo
         x(i) = sum/a(p(i),i)
    enddo
         solver = x


endfunction
!-------------------------------------------------------------
subroutine update
implicit none
real(8) :: rhov,rhol,pi,ci,yl,yv,rhoi,sint
real(8) :: drvdp,drldp,drdp,drdyv,dv(3)
real(8) :: a(3,3),b(3)
res = 0.0
!compute the jacobian solve and update
 do i=1,nx
        
    sint = 0.5d0*(S(i) + S(i-1))
    rhov = max(v(2,i),vp0)/(Rvap*t0)
    rhol = rhow(v(2,i))
      yv = v(1,i)
      yl = 1.0 - yv 
      pi = v(2,i)
    rhoi = rho(i)
    
   drvdp = 1.0/(Rvap*t0)           
   drldp = 1.0/(k0*(1.0d0 +(n0/k0)*(pi-p0))/(rhol))
     
    drdp = (rhoi**2)*(yv/((rhov**2)*(1.0/drvdp)) + yl/((rhol**2)*(1.0/drldp)))

   drdyv = -1.0*(rhoi**2)*(1.0d0/rhov - 1.0d0/rhol)
            
 dudv(:,:,i) = 0.0d0  

 dudv(1,1,i) = rhoi  + yv*drdyv
 dudv(1,2,i) = drdp*yv
 dudv(1,3,i) = 0.0d0 
 
 dudv(2,1,i) = drdyv
 dudv(2,2,i) = drdp
 dudv(2,3,i) = 0.0d0

 dudv(3,1,i) = v(3,i)*drdyv
 dudv(3,2,i) = drdp*v(3,i)
 dudv(3,3,i) = rhoi

      rhs(:,i) = sint*Sm(:,i) - (S(i)*F(:,i) - S(i-1)*F(:,i-1))/(h)
      res = res + rhs(:,i)**2
      a = sint/dt*dudv(:,:,i)
      b = rhs(:,i)
      !dv = 0.0d0
      dv = solver(a,b)
     
      !print*, dudv(1,:,i)      
      !print*, dudv(2,:,i)      !stop
      !print*, dudv(3,:,i)    
      !print*,dv(:)
      !print*,rhs(:,i),Sm(:,i)
      !print*, sint
     v(:,i) = v(:,i) +  dv 
 enddo
 !stop
     !do i = 0,nx
        !print*, F(:,i)
     !enddo
     !stop
     res = sqrt(res)
endsubroutine
!-----------------------------------------------------------
subroutine boundary
implicit none
real(8) :: mi,yvi,pi,ui,rhoi,ci
real(8) :: yv,yl,rhov,rhol
integer :: j
! inlet condition
   v(1,0) = yv0  
   v(3,0) = v(3,1)
       
      yv = v(1,0)
      
do j=1,1000                                             
      v(2,0) = P0stag - 0.5d0*rho(0)*(v(3,0)**2)
      rhol = rhow(v(2,0)) 
      rhov = max(vp0,v(2,0))/(Rvap*t0)
      yl = 1.0 - yv
     rho(0)= 1.0/(yv/rhov+yl/rhol)  
enddo                                                   
   c(0)   = sound(v(2,0),v(1,0),rho(0)) 


!print*, v(2,1)       
!outlet condition

   mi = v(3,nx)/c(nx)

  if(mi<1.0) then
    v(2,nx+1) = Pback
 else if (mi>1.0) then
    v(2,nx+1) = v(2,nx)
 else 
    print*, 'Incorrect Mach Number',mi,iter
    stop
 endif

    v(1,nx+1) = v(1,nx)
    v(3,nx+1) = v(3,nx)


endsubroutine
!-----------------------------------------------------------
subroutine source
implicit none
integer :: i
real(8) :: p,rhov,rhol,yl,yv,dp,eq,sint
real(8),parameter :: ctani = 1.0d0
real(8),parameter :: ktani = 5.0d0

    Sm = 0.0d0

do i=1,nx
   p = v(2,i)
 sint= (S(i)+S(i-1))*0.5d0
   if(p<vp0) then
!    print*, 'hi' 
     rhov = max(p,vp0)/(Rvap*t0) 
     rhol = rhow(p)

      dp = min(vp0-p,ktani*vp0)
      eq = sqrt(2.0/3.0*dp/rhol)      

      yv = v(1,i)

      Sm(1,i) = ctani*dp*(1.0-yv)*rhov*eq
         !print*, Sm(:,i)
   endif
      Sm(3,i) = v(2,i)*(S(i)-S(i-1))/h/sint
enddo
 
endsubroutine
!-----------------------------------------------------------
subroutine initial
implicit none
integer :: i,j,n
real(8) :: dr,yl,yv,rhol,rhov


rho = r0
 yv =yv0
 
v(1,:) = yv
v(2,:) = P0stag
 do i=0,nx+1
    v(3,i) = u0/(S(i)/2.0) 
 enddo
do j=1,1000                                              
   dr =0.0d0
   do i=0,nx+1        
      v(2,i) = P0stag - 0.5d0*rho(i)*(v(3,i)**2)
          dr = dr + dabs(rhow(v(2,i))-rho(i))
      rhol = rhow(v(2,i)) 
      rhov = max(vp0,v(2,i))/(Rvap*t0)
      yl = 1.0 - yv
     rho(i)= 1.0/(yv/rhov+yl/rhol)  
!       print*, j          
   enddo
   
   if(dr.le.1d-4) exit
enddo
  pin = v(2,0)
!print*, S
!stop
v(1,nx+1) = yv0 
!v(2,nx+1) = Pback
!v(3,nx+1) = v(3,nx)
 
endsubroutine
!-----------------------------------------------------------
subroutine therm
implicit none
real(8) :: rhov,rhol,yv,yl,rhoi
real(8) :: pi,eq
integer :: i

do i = 0,nx+1

     yv = v(1,i)
     pi = v(2,i)
   rhov = max(pi,vp0)/(Rvap*t0)
   rhol = rhow(pi)
     
 yl = 1.0 - yv   
      eq = yv/rhov + yl/rhol
   rhoi = 1.0/(eq)  
   rho(i) = rhoi
     c(i) = sound(pi,yv,rhoi)
  !print*, pi
enddo
!stop
endsubroutine
!-----------------------------------------------------------
subroutine flux
implicit none
integer :: i,j
real(8) :: CVLp,CVLn,Dp,Dn,Dpi,Dni1,api,bi,bi1,ani1,machpi,machni1,machi,machi1, &
           Fci(3),Fci1(3),Pi(3),Pi1(3),CLSp,CLSn,CEp,CEn,M_1,Mp,Mn,delp
real(8) :: Pvmod(3),Pemod(3),c1
real(8) :: cint,ui,ui1,p_i,p_i1,rhoi,rhoi1 

  Pemod =(/0.0d0,0.0d0,1.0d0/)

 if(iflux ==1) then  !unmodified Luoi stefen
   do i=0,nx
      
    cint = sqrt((rho(i)*(c(i)**2) + rho(i+1)*(c(i+1)**2))/(rho(i)+rho(i+1)))
     ui  = v(3,i)
     ui1 = v(3,i+1)
  machi  = ui/(cint)
  machi1 = ui1/(cint)
    rhoi = rho(i)
    rhoi1= rho(i+1)
      p_i = v(2,i)
      p_i1= v(2,i+1)

    api= 0.5d0*(1.0d0+sign(1.0d0,machi))
    ani1= 0.5d0*(1.0d0-sign(1.0d0,machi1))
     bi= -max(0.0d0,1.0d0-int(dabs(machi)))
    bi1= -max(0.0d0,1.0d0-int(dabs(machi1)))

    machpi= 0.25d0*(machi + 1.0d0)**2
   machni1=-0.25d0*(machi1- 1.0d0)**2
       Dpi=0.25d0*((machi+1.0)**2)*(2.0-machi)
      Dni1=0.25d0*((machi1-1.0)**2)*(2.0+machi1)

    CVLp = api*(1.0d0+bi)*machi - bi*machpi
    CVLn = ani1*(1.0d0+bi1)*machi1-bi1*machni1

    Dp   = api* (1.0+bi)-bi*Dpi
    Dn   = ani1*(1.0+bi1)-bi1*Dni1

    Fci = (/ v(1,i),1.0d0,ui    /) 
    Fci1= (/ v(1,i+1),1.0d0,ui1 /)

    Pi  =(/0.0d0,0.0d0,p_i/)
    Pi1 =(/0.0d0,0.0d0,p_i1/)       

    CLSp =max(0.0,CVLp+CVLn)
    CLSn =min(0.0,CVLp+CVLn)
 
  F(:,i)=(cint)*(rhoi*CLSp*Fci + rhoi1*CLSn*Fci1) + Dp*Pi + Dn*Pi1
  

    enddo     
!stop
 elseif (iflux == 2) then !unmodified LDFSS
   do i=0,nx

    cint = sqrt((rho(i)*(c(i)**2) + rho(i+1)*(c(i+1)**2))/(rho(i)+rho(i+1)))
     ui  = v(3,i)
     ui1 = v(3,i+1)
  machi  = ui/(cint)
  machi1 = ui1/(cint)
    rhoi = rho(i)
    rhoi1= rho(i+1)
      p_i = v(2,i)
      p_i1= v(2,i+1)

    api= 0.5d0*(1.0d0+sign(1.0d0,machi))
    ani1= 0.5d0*(1.0d0-sign(1.0d0,machi1))
     bi= -max(0.0d0,1.0d0-int(dabs(machi)))
    bi1= -max(0.0d0,1.0d0-int(dabs(machi1)))

    machpi= 0.25d0*(machi + 1.0d0)**2
   machni1=-0.25d0*(machi1- 1.0d0)**2
       Dpi=0.25d0*((machi+1.0)**2)*(2.0-machi)
      Dni1=0.25d0*((machi1-1.0)**2)*(2.0+machi1)

    CVLp = api*(1.0d0+bi)*machi - bi*machpi
    CVLn = ani1*(1.0d0+bi1)*machi1-bi1*machni1

    Dp   = api* (1.0+bi)-bi*Dpi
    Dn   = ani1*(1.0+bi1)-bi1*Dni1

    Fci = (/ v(1,i),1.0d0,ui    /) 
    Fci1= (/ v(1,i+1),1.0d0,ui1 /)

    Pi  =(/0.0d0,0.0d0,p_i/)
    Pi1 =(/0.0d0,0.0d0,p_i1/)                                          
    M_1 = 0.25*bi*bi1*(sqrt(0.5*(machi**2+machi1**2))-1.0d0)**2
    delp= p_i-p_i1

     Mp= M_1*max(0.0,1.0-((delp/(p_i+p_i1))+(2.0*delp/p_i)))
     Mn= M_1*max(0.0,1.0+((delp/(p_i+p_i1))-(2.0*delp/p_i1)))

     CEp= CVLp - Mp
     CEn= CVLn + Mn

  F(:,i)=(cint)*(rhoi*CEp*Fci + rhoi1*CEn*Fci1) + Dp*Pi + Dn*Pi1

   !print*, Mn,CEn
   enddo 
 elseif (iflux == 3) then !modified Loiu Stefan
   do i=0,nx

    cint = sqrt((rho(i)*(c(i)**2) + rho(i+1)*(c(i+1)**2))/(rho(i)+rho(i+1)))
     ui  = v(3,i)
     ui1 = v(3,i+1)
  machi  = ui/(cint)
  machi1 = ui1/(cint)
    rhoi = rho(i)
    rhoi1= rho(i+1)
      p_i = v(2,i)
      p_i1= v(2,i+1)
   !print*, cint                                                                    
    api= 0.5d0*(1.0d0+sign(1.0d0,machi))
    ani1= 0.5d0*(1.0d0-sign(1.0d0,machi1))
     bi= -max(0.0d0,1.0d0-int(dabs(machi)))
    bi1= -max(0.0d0,1.0d0-int(dabs(machi1)))
                                                                      
    machpi= 0.25d0*(machi + 1.0d0)**2
   machni1=-0.25d0*(machi1- 1.0d0)**2
       Dpi=0.25d0*((machi+1.0)**2)*(2.0-machi)
      Dni1=0.25d0*((machi1-1.0)**2)*(2.0+machi1)
                                                                      
    CVLp = api*(1.0d0+bi)*machi - bi*machpi
    CVLn = ani1*(1.0d0+bi1)*machi1-bi1*machni1
                                                                      
    Dp   = api* (1.0+bi)-bi*Dpi
    Dn   = ani1*(1.0+bi1)-bi1*Dni1
                                                                      
    Fci = (/ v(1,i),1.0d0,ui    /) 
    Fci1= (/ v(1,i+1),1.0d0,ui1 /)
                                                                      
    Pi  =(/0.0d0,0.0d0,p_i/)
    Pi1 =(/0.0d0,0.0d0,p_i1/)                                         

                  
    Pvmod = 0.5d0*(Pi+Pi1) +0.5d0*(Dp-Dn)*(Pi-Pi1) +0.5d0*Pemod*(rhoi+rhoi1)*(Dp+Dn-1.0d0)*(cint**2)
 
   CLSp =max(0.0,CVLp+CVLn)
   CLSn =min(0.0,CVLp+CVLn)  

  F(:,i)=(cint)*(rhoi*CLSp*Fci + rhoi1*CLSn*Fci1) + Pvmod

  !print*, CLSP,CLSn    !print*, CLSP,CLSn
  ! print*, i,F(:,i)         
   enddo 
 elseif (iflux == 4) then !modified LDFSS
   do i=0,nx
 
     cint =sqrt((rho(i)*(c(i)**2) + rho(i+1)*(c(i+1)**2))/(rho(i)+rho(i+1)))
      ui  = v(3,i)
      ui1 = v(3,i+1)
   machi  = ui/(cint)
   machi1 = ui1/(cint)
     rhoi = rho(i)
     rhoi1= rho(i+1)
       p_i = v(2,i)
       p_i1= v(2,i+1)
 
     api= 0.5d0*(1.0d0+sign(1.0d0,machi))
     ani1= 0.5d0*(1.0d0-sign(1.0d0,machi1))
      bi= -max(0.0d0,1.0d0-int(dabs(machi)))
     bi1= -max(0.0d0,1.0d0-int(dabs(machi1)))
 
     machpi= 0.25d0*(machi + 1.0d0)**2
    machni1=-0.25d0*(machi1- 1.0d0)**2
        Dpi=0.25d0*((machi+1.0)**2)*(2.0-machi)
       Dni1=0.25d0*((machi1-1.0)**2)*(2.0+machi1)
 
     CVLp = api*(1.0d0+bi)*machi - bi*machpi
     CVLn = ani1*(1.0d0+bi1)*machi1-bi1*machni1
 
     Dp   = api* (1.0+bi)-bi*Dpi
     Dn   = ani1*(1.0+bi1)-bi1*Dni1
 
     Fci = (/ v(1,i),1.0d0,ui    /) 
     Fci1= (/ v(1,i+1),1.0d0,ui1 /)
 
     Pi  =(/0.0d0,0.0d0,p_i/)
     Pi1 =(/0.0d0,0.0d0,p_i1/)                                         
     M_1 = 0.25*bi*bi1*(sqrt(0.5*(machi**2+machi1**2))-1.0d0)**2
     delp= p_i-p_i1
 
       Mp= M_1*max(0.0,1.0-0.5d0*((delp+dabs(delp))/(rhoi*(cint)**2)))
       Mn= M_1*max(0.0,1.0+0.5d0*((delp-dabs(delp))/(rhoi1*(cint)**2)))
  
       CEp= CVLp - Mp
       CEn= CVLn + Mn
  Pvmod = 0.5d0*(Pi+Pi1) +0.5d0*(Dp-Dn)*(Pi-Pi1) +0.5d0*Pemod*(rhoi+rhoi1)*(Dp+Dn-1.0d0)*(cint**2)
      F(:,i)=(cint)*(rhoi*CEp*Fci + rhoi1*CEn*Fci1) + Pvmod
      !print*, Dp*Pi+Dn*Pi1
   !print*, i,F(:,i)
   !print*, Cep,Cen
   enddo
   !stop
 else
    print*, 'Incorrect flux flag selected'
 endif


endsubroutine
!-------------------------------------------------------------------------------
!Tait equation
function rhow(press)    
 real(8) :: press,rhow
 real(8) :: a
 a = (1.0+n0/k0*(press-p0))  
 rhow = r0*(a)**(1.0/n0)
 return
 endfunction rhow
!---------------------------------------------------------------------------------

function sound(p,yv,rhoi)   
real(8) :: p,yv,rhoi,sound,yl
real(8) :: cv,cl,av,al,rhol,rhov
   
 rhol = rhow(p)
 rhov = max(p,vp0)/(Rvap*t0)
 
   yl = 1.0 - yv
   cv = sqrt(Rvap*t0)            ! 1/drhovdp
   cl = sqrt(k0*(1.0d0 +(n0/k0)*(p-p0))/(rhol))

   av = (rhoi-rhol)/(rhov-rhol)
   al = 1.0d0 - av
 
 sound = rhov*(cv**2)*rhol*(cl**2)/(rhoi*(av*rhol*(cl**2)+al*rhov*(cv**2)))     
 
 sound = sqrt(sound)   
 !print*, cv,cl,p,rhol
endfunction sound

end program 
