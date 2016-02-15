!------------------------------------------------------------------------------------
subroutine PBCsmear01(drsq, a_pol12, aewald, ts0, ts1)
use system_mod
use math
implicit none
double precision, parameter :: g23=1.3541179394264d0
double precision, parameter :: TSP= 1.12837916709551257390d0 !  2/sqrt(pi)

double precision :: dd,dri,drsqi,expon,AA, a, drsq, aewald, rA,rA3,exp1,a_sqrt3, ch1,ch2,sr1,sr2,Rc_tholesq, ts0, ts1, a_pol12
double precision :: er
!double precision, external :: gammq2_3 !changed by D. Elton. These are now in the 'math' module
!double precision, external :: erfc

dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri
expon = dexp(-(aewald*aewald*dRsq))
   
er = erfc(aewald*dd)
ts0=er*dRi
ts1 = (ts0 + const_ts1*expon)*drsqi

Rc_tholesq=25.d0
if (dRsq<Rc_tholesq) then
   rA3 = dd*drsq*a_pol12
   ra = a_pol12**(1.d0/3.d0)
   exp1 = dexp(-rA3)
   ch1 = -1.d0*dri*drsqi
   sr1 = (1.d0 -exp1)*dri*drsqi
   ts1 = ts1+ch1+sr1
   ts0 = ts0  - exp1*dri + rA*g23*gammq2_3( ra3 )
endif

end subroutine PBCsmear01
!------------------------------------------------------------------------------------
subroutine PBCsmear2(drsq, a_pol12, aewald, ts1, ts2)
use system_mod
use consts
use math
implicit none
double precision, intent(in) :: drsq, a_pol12, aewald
double precision, intent(out) :: ts1, ts2
double precision :: dd,dri,drsqi,expon,AA, rA,rA3,exp1,ch1,ch2,sr1,sr2, Rc_tholesq
double precision :: er
!double precision, external :: erfc


dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri
expon = dexp(-(aewald*aewald*dRsq))

er = erfc(aewald*dd)
ts1 = (er*dri + const_ts1*expon)*drsqi
ts2 = (ts1    + const_ts2*expon)*drsqi

Rc_tholesq=25.d0
if (dRsq<Rc_tholesq) then
   rA3 = (dd*drsq)*a_pol12!**3
   exp1 = dexp(-rA3)
   ch1 = -1.d0*dri*drsqi
   ch2 = ch1*drsqi
   sr1 = (1.d0 -exp1)*dri*drsqi
   sr2 = (sr1 - exp1*a_pol12) *drsqi
   ts1 = ts1+ch1+sr1
   ts2 = ts2+ch2+sr2
endif

end subroutine PBCsmear2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PBCsmear3(drsq, a_pol12,  aewald, ts1, ts2, ts3)
use system_mod
use consts
use math
implicit none
double precision, intent(in) :: drsq, a_pol12,  aewald
double precision, intent(out) :: ts1, ts2, ts3
double precision :: dd,dri,drsqi,expon,AA,rA,rA3,exp1,ch1,ch2,ch3,sr1,sr2,sr3,Rc_tholesq
double precision :: er
!double precision, external :: erfc

dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri
expon = dexp(-(aewald*aewald*dRsq))

er = erfc(aewald*dd)
ts1 = (er*dri + const_ts1*expon)*drsqi
ts2 = (ts1    + const_ts2*expon)*drsqi
ts3 = (ts2    + const_ts3*expon)*drsqi

Rc_tholesq=25.d0
if (dRsq<Rc_tholesq) then
   rA3 = (dd*drsq)*a_pol12!**3
   exp1 = dexp(-rA3)
   ch1 = -1.d0*dri*drsqi
   ch2 = ch1*drsqi
   ch3 = ch2*drsqi
   sr1 = (1.d0 -exp1)*dri*drsqi
   sr2 = (sr1 - exp1*a_pol12) *drsqi
   sr3 = (sr2 - 0.6d0*exp1*dd*a_pol12*a_pol12) * drsqi
   ts1 = ts1+ch1+sr1
   ts2 = ts2+ch2+sr2
   ts3 = ts3+ch3+sr3
endif

end subroutine PBCsmear3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine self1(drsq, aewald, ts0, ts1)
use system_mod
use consts
use math
implicit none
double precision, intent(in) :: drsq
double precision, intent(out) :: ts0, ts1
double precision :: dd, dri, expon, drsqi, kr0, kr1, ch0, ch1, aewald
!double precision, external :: erfc

dd = dsqrt(drsq)
drsqi = 1.d0/drsq
dri = dsqrt(drsqi)
expon = dexp(-(aewald*aewald*drsq))
kr0=(erfc(aewald*dd))*dRi
kr1=(kr0 + const_ts1*expon)*drsqi
ch0 = -dri
ch1 = ch0*drsqi
ts0 = ch0 + kr0
ts1 = kr1 + ch1

end subroutine self1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine self2(drsq, aewald, ts1, ts2)
use system_mod
use consts
use math
implicit none
double precision, intent(in) :: drsq, aewald
double precision, intent(out) :: ts1, ts2
double precision :: dd, dri, expon, drsqi, kr0, kr1, kr2, ch0, ch1, ch2
!double precision, external :: erfc

dd = dsqrt(drsq) 
drsqi = 1.d0/drsq 
dri = dsqrt(drsqi)
expon = dexp(-(aewald*aewald*drsq))
kr0=(erfc(aewald*dd))*dRi
kr1=(kr0 + const_ts1*expon)*drsqi
kr2=(kr1 + const_ts2*expon)*drsqi
ch0 = -dri
ch1 = ch0*drsqi
ch2 = ch1*drsqi
ts1 = kr1 + ch1
ts2 = kr2 + ch2

end subroutine self2

