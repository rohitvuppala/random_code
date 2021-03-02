
!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!
!!!! O2- species 1;O - species 2;!!!!!!!!!!!!
!!!!RK-2 Method!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program main
implicit none

real(8),parameter :: Runiv=8314.34d0,M1=32d0,M2=16d0,TB1=1d0,TB2=3d0,Cf=2.75e16,eta=-1d0, &
                     theta=59500d0,Ti= 4150,rhoi=3.35e-4                                                              !intial,fixed values

real(8),parameter,dimension(7) :: A1=(/3.660960830d+00,6.56365523d-04,-1.41149485d-07, &                   ! O2 aboveT=1000K
                                       2.05797658d-11,-1.29913248d-15,-1.21597725d+03,3.41536184d+00/), &
                                  A2=(/2.54363697d+00,-2.73162486d-05,-4.19029520d-09, &                   ! O  aboveT=1000K
                                       4.95481845d-12,-4.79553694d-16,2.92260120d+04,4.92229457d+00/)
real(8),parameter,dimension(7) :: B1=(/3.78245636d+00,-2.99673415d-03,9.84730200d-06,  &
                                      -9.68129508d-09,3.24372836d-12,-1.06394356d+03,3.65767573d+00/), &    ! O2 below T=1000K
                                  B2=(/3.16826710d+00,-3.27931884d-03,6.64306396d-06,  &
                                      -6.12806624d-09,2.11265971d-12,2.91222592d+04,2.05193346d+00/)        ! O below T=1000K
real(8),dimension(7):: A,B,A1h,B1h,A2h,B2h

real(8) :: T,rho1,rho2,kf,hi,ei,e1,e2,engi,eng,eng1,eng2,h1,h2,d1dt,gp,gr,dg,keq,dt,deng,dengdT,time,k1,k2, &
           rho1sp,rho2sp,res,res0

integer :: nitermax,nmax,n,i
rho1=0d0
rho2=0d0
T=Ti
dt= 1e-6
nitermax = 1000
nmax = 100000000
time=0

! Cp/Runiv formula for curve fitting
! Cp/Runiv = a(1)+a(2)*T+a(3)*T**2+a(4)*T**3+a(5)*T**4+a(6)+a(7)
!!!! a = [A,B,C,D,E,F,G], for Enthalpy/Runiv calculation, divide by molecular weight for the constants

 A1h=A1/32;
 A2h=A2/16;
 B1h=B1/32;
 B2h=B2/16;

!!!!!!! intial condition
!!!h/R(molar)=Cp/Runiv = a(1)*T+a(2)*T**2/2+a(3)*T**3/3+a(4)*T**4/4+a(5)*T**5/5+a(6)
!!!h/R(mass) we divided coefficients  above by Molecular weight
rho1=rhoi
rho2=rhoi-rho1
hi =Runiv*(A1h(1)*Ti        &
      +A1h(2)*(Ti**2)/2     &
      +A1h(3)*(Ti**3)/3     &
      +A1h(4)*(Ti**4)/4     &
      +A1h(5)*(Ti**5)/5     &
      +A1h(6))

ei  = hi - Runiv*Ti/M1
engi=rhoi*ei

!print*, engi

!!!!!!!!!!!
!!!!!!!check if the temperature is above or below 1000K!!!!!!
!!!!!!!!!!!

open (unit = 2, file = "data.dat")
!!!Loop starts
!!!!!!!!!
!!!!!!!!!
do i=1,nmax

if (T<=1000) then
A=B1h
B=B2h
elseif (T>1000) then
A=A1h
B=A2h
endif

!!!!!!keq,kf computation

gr = Runiv*M1*(A(1)*(1-log(T))*T         &
   -A(7)*T                        &
   -A(2)*(T**2)/2                 &
   -A(3)*(T**3)/6                 &
   -A(4)*(T**4)/12                &
   -A(5)*(T**5)/20+A(6))


gp = Runiv*2*M2*(B(1)*(1-log(T))*T      &
   -B(7)*T                        &
   -B(2)*(T**2)/2                 &
   -B(3)*(T**3)/6                 &
   -B(4)*(T**4)/12                &
   -B(5)*(T**5)/20+B(6))

dg = gp-gr



keq =(101325d0/Runiv/T)*exp(-dg/Runiv/T)
kf = (2.75d16/T)* exp(-59500d0/T)


d1dt= (-1)*M1*(kf)*((rho1/M1)- rho2*rho2/(M2*M2*keq))*(TB1*rho1/M1 + TB2*rho2/M2)

rho1sp = rho1 + dt*d1dt
rho2sp = rhoi - rho1
d1dt= (-1)*M1*(kf)*((rho1sp/M1)- rho2sp*rho2sp/(M2*M2*keq))*(TB1*rho1sp/M1 + TB2*rho2sp/M2)
rho1 = 0.5* rho1sp + 0.5*(rho1 + d1dt*dt )


rho2=rhoi-rho1


!!!!!!!!!!!!!
!!Net energy change = 0
!!!!!!!!!!!!!!

do n = 1,nitermax

if (T<=1000) then
  A=B1h
  B=B2h
elseif (T>1000) then
  A=A1h
  B=A2h
endif

h1 =  (A(1)*T+A(2)*(T**2)/2 + A(3)*(T**3)/3+A(4)*(T**4)/4+ &
       A(5)*(T**5)/5+A(6))*Runiv

e1 = h1- Runiv*T/M1
eng1=rho1*e1

h2 =  (B(1)*T+B(2)*(T**2)/2 + B(3)*(T**3)/3+B(4)*(T**4)/4+ &
       B(5)*(T**5)/5+B(6))*Runiv

e2 = h2- Runiv*T/M2
eng2=rho2*e2

eng=eng1+eng2
deng=eng-engi


!!!!!!!!! Update T to make deng =0

dengdT = (rho1*( (A(1)+A(2)*(T) + A(3)*(T**2)+A(4)*(T**3)+                     &
       A(5)*(T**4)) - 1/M1) + rho2*(B(1)+B(2)*(T)/2 + B(3)*(T**2)+B(4)*(T**3)+ &
       B(5)*(T**4) - 1/M2))*Runiv

T = T - deng/dengdT


if (abs(deng/dengdT) < 10e-6) then

   exit

end if

end do

print*, i,n,time,T,rho1/rhoi,dg

write(2,*) i,n,time,rho1/rhoi,rho2/rhoi,T
if (i.eq.1) res0=abs(d1dt)
res = abs(d1dt)/res0
if ((abs(res)<10e-8)) then
   exit
end if
time = time + dt
end do

close(2)

end










