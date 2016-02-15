subroutine pot_spc(RR, En, dRR, virt, dip_momI, chg)
use consts
use system_mod
use pot_mod
use neigh_mod
use math
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR
double precision, intent(out) :: En 
double precision, dimension(3, Natoms), intent(out) :: dRR
double precision, dimension(3, 3), intent(out) :: virt
double precision, dimension(3, NWaters), intent(out) ::  dip_momI
double precision, dimension(Natoms), intent(out) ::  chg
double precision, dimension(3) :: dip_mom
integer :: itmp, i, is, js, iw, jw, iO, ih1, ih2,  iOa, iH1a, iH2a,  jO, jH1, jH2
integer :: iat, jat, isp, jsp, jsp0, i3, j3, ii, jj, kk, ix, iy, j
double precision :: Ureck, e1, fac
double precision :: tmp, R2i, R4i, R6i, R8i, uij, dRsq, dRij, b2
double precision :: dri, drsqi, expon, er, dd, ch3, exp1, qi, qj
!double precision :: Umon, Uvdw, Uelec
double precision, dimension(3) ::  Ri,Rij, roh1,roh2, rij0, dh1, dh2, Rhh,rh1O,rh2O
double precision, dimension(3, 3) :: r1, dr1, vij
double precision, dimension(3), save :: pr_box
integer , save :: pr_Natoms
double precision :: dROH1, dROH2, dAxB, costh, sinth, ang, tmp1, tmp2, dvdd1, dvdd2, dvdd3, rhh0, drhh, dd1, dd2, dd3
double precision ::  kr0, kr1, ch0, ch1, ts0, ts1
double precision, dimension(3,3) :: Atmp
double precision, dimension(3) :: AxB, dij, tmp_der
double precision, dimension(5) :: tmp_arr5
!double precision, external :: erfc
logical :: remove_net_force
!type t_neigh 
!   integer :: N
!   integer, dimension(2, 2000) :: ij
!   double precision, dimension(3, 2000) :: Rij
!   double precision, dimension(2000) :: R2
!end type t_neigh 
type(t_neigh) :: neigh
!!!!!!!!!!!!!!!!!!!!!!!!!!
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.not. allocated(R) .or. Natoms/=pr_Natoms .or.  &
      dsqrt(dot_product(box-pr_box, box-pr_box))>1.d-12) then
   pr_Natoms = Natoms
   pr_box = box
endif

virt = 0.d0
vir_rec = 0.d0
Umon = 0.d0; Uvdw = 0.d0; Uelec = 0.d0
phi= 0.d0
Efq = 0.d0
dR = 0.d0
do iw=1, Nwaters
   iOa=3*iw-2; iH1a=3*iw-1; iH2a=3*iw
   iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters;
   R(1:3, (/iO, iH1, iH2/)) = RR(1:3, (/iOa, ih1a, ih2a/) )! * boxi
   chg( (/ioa, ih1a, ih2a/) ) = charge( (/io, ih1, ih2/) )
enddo

!call find_neigh(1, R, neigh)
!print*,'N = ',neigh % N
!do j=1, neigh % N
!   write(*,'(i4,2x,i3, 4(f12.6,1x))')j, neigh % j(j), neigh % Rij(1:3, j), dsqrt(neigh % R2(j))
!enddo


if (pot_model==4) then !  Harmonic oscillators (qspcfw)
 do iw=1, Nwaters
      iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters
      r1(1:3, 1:3) = R(1:3, (/iO, ih1, ih2/) )
!...............
      ROH1(:)=r1(1:3,2)-r1(1:3,1);roh1=roh1-box*anint(roh1/box);dROH1=dsqrt(dot_product(ROH1, ROH1))
      ROH2(:)=r1(1:3,3)-r1(1:3,1);roh2=roh2-box*anint(roh2/box);dROH2=dsqrt(dot_product(ROH2, ROH2))
      AxB(1) = ROH1(2)*ROH2(3) - ROH1(3)*ROH2(2)
      AxB(2) =-ROH1(1)*ROH2(3) + ROH1(3)*ROH2(1)
      AxB(3) = ROH1(1)*ROH2(2) - ROH1(2)*ROH2(1)
      dAxB = dsqrt(dot_product(AxB, AxB) )
      costh = dot_product(ROH1, ROH2) / (dROH1*dROH2)
      sinth = dAxB / (dROH1*dROH2)
      ang = atan2(sinth, costh)
      e1 = harmkb*( (dROH1-roh0)**2 + (dROH2-roh0)**2) + harmka*(ang-th0)**2
      tmp1 = 2.d0*harmkb*(dROH1-roh0) / dROH1
      tmp2 = 2.d0*harmkb*(dROH2-roh0) / dROH2
      tmp = 2.d0*harmka*(ang-th0)
      Atmp(1:3,1)= -(ROH2(1:3)/(dROH1*dROH2) - costh*ROH1(1:3)/(dROH1*dROH1))/sinth
      Atmp(1:3,2)= -(ROH1(1:3)/(dROH1*dROH2) - costh*ROH2(1:3)/(dROH2*dROH2))/sinth
      Atmp(1:3,3) = -Atmp(1:3,1)-Atmp(1:3,2)
      dh1(1:3) = tmp1*ROH1(:) + tmp*Atmp(:,1)
      dh2(1:3) = tmp2*ROH2(:) + tmp*Atmp(:,2)
      dR1(:,2) = dh1(1:3)
      dR1(:,3) = dh2(1:3) 
      dR1(:,1) =-dh1(1:3) - dh2(1:3)
      Umon = Umon + e1
      dR(1:3, (/iO, ih1, ih2/)) = dr1

      do iy=1,3
         virt(1:3,iy)=virt(1:3,iy) + ROH1(1:3)*dh1(iy) + ROH2(1:3)*dh2(iy)
      enddo


     ! virialc = virialc + dot_product(ROH1,dh1)
      !virialc = virialc + dot_product(ROH2,dh2)

   enddo
else if (pot_model==5) then   ! More oscillators (spcfw) JCP 106 2400
 rhh0 = 2.d0*dsin(0.5d0*th0)*roh0
 do iw=1, Nwaters
      iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters
      r1(1:3, 1:3) = R(1:3, (/iO, ih1, ih2/) )
!...............
      ROH1(:)=r1(1:3,1)-r1(1:3,2);roh1=roh1-box*anint(roh1/box);dROH1=dsqrt(dot_product(ROH1, ROH1))!PBC
      ROH2(:)=r1(1:3,1)-r1(1:3,3);roh2=roh2-box*anint(roh2/box);dROH2=dsqrt(dot_product(ROH2, ROH2))!PBC
      RHH(:) =r1(1:3,2)-r1(1:3,3);rhh =rhh -box*anint(rhh /box);dRhh =dsqrt(dot_product(Rhh , Rhh ))!PBC
      dROH1=dsqrt(ROH1(1)*ROH1(1) + ROH1(2)*ROH1(2) + ROH1(3)*ROH1(3) )
      dROH2=dsqrt(ROH2(1)*ROH2(1) + ROH2(2)*ROH2(2) + ROH2(3)*ROH2(3) )
      dRHH =dsqrt(RHH (1)*RHH (1) + RHH (2)*RHH (2) + RHH (3)*RHH (3) )
      dd1 = dROH1 - roh0
      dd2 = dROH2 - roh0
      dd3 = dRHH  - rhh0
      e1 = spcf_Ap*(dd1*dd1 + dd2*dd2) + spcf_Bp*dd3*dd3 + spcf_Cp*(dd1+dd2)*dd3 + spcf_Dp*dd1*dd2

      dVdd1 =2.d0*spcf_Ap*dd1 + spcf_Cp*dd3 + spcf_Dp*dd2
      dVdd2 =2.d0*spcf_Ap*dd2 + spcf_Cp*dd3 + spcf_Dp*dd1
      dVdd3 =2.d0*spcf_Bp*dd3 + spcf_Cp*(dd1+dd2)

      dR1(1:3,1) =  dVdd1*ROH1/dROH1 + dVdd2*ROH2/dROH2
      dR1(1:3,2) = -dVdd1*ROH1/dROH1 + dVdd3*RHH/dRHH
      dR1(1:3,3) = -dVdd2*ROH2/dROH2 - dVdd3*RHH/dRHH
      Umon = Umon + e1
      dR(1:3, (/iO, ih1, ih2/)) = dr1
      do iy=1,3
!         virt(1:3,iy)=virt(1:3,iy) + ROH1(1:3)*dR1(iy) + ROH2(1:3)*dR2(iy)
         virt(1:3,iy)=virt(1:3,iy) - ROH1(1:3)*dR1(iy,2) - ROH2(1:3)*dR1(iy,3)
      enddo
   enddo
endif
!*
!* Calculate the pairwise additive VDW and Electrostatic (Real) Interactions
!*
   Uelec = 0.d0
   do i=1, Nwaters
      iO = i
      call find_neigh(iO, R, Neigh)
      do j=1, neigh % N
         jO = neigh % j(j)
         Rij = neigh % Rij(1:3, j)
         dRsq = neigh % R2(j)
      !....... VDW ENERGY........!
         dRij=dsqrt(dRsq);
         R2i = 1.d0/dRsq; R4i=R2i*R2i; R6i=R2i*R2i*R2i; R8i=R6i*R2i
    
         uij = r6i*(vdwA*r6i+vdwC)
         dij = -R8i*(12.d0*vdwA*R6i + 6.d0*vdwC)*Rij
         Uvdw = Uvdw + uij
         dR(1:3,iO) = dR(1:3,iO) + dij
         dR(1:3,jO) = dR(1:3,jO) - dij
         virt(1:3,1) = virt(1:3,1) + dij(1)*Rij
         virt(1:3,2) = virt(1:3,2) + dij(2)*Rij
         virt(1:3,3) = virt(1:3,3) + dij(3)*Rij
         !...... Electrostaric (REAL)........!
         Rij0 = Rij - R(1:3, iO) + R(1:3, jO) !???
         do isp=1, 3
            iat = iO +  (isp-1)*Nwaters
            Ri = R(1:3, iat) + rij0
            qi = charge(iat)
            do jsp=1, 3
               jat = jO +  (jsp-1)*Nwaters
               qj = charge(jat)
               Rij = Ri - R(1:3, jat) 
	       Rij = Rij - box*anint(Rij*boxi) !PBC
               drsq = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
               dd  = dsqrt(drsq)
               expon = dexp(-(aewald*aewald*dRsq))
               ts0 = erfc(aewald*dd) / dd
               ts1 = (ts0 + const_ts1*expon)/drsq
    
               phi(iat) = phi(iat) + ts0*qj
               phi(jat) = phi(jat) + ts0*qi
               Uelec = Uelec + ts0*qi*qj
               Efq(1:3,iat) = Efq(1:3,iat) + ts1*qj*Rij
               Efq(1:3,jat) = Efq(1:3,jat) - ts1*qi*Rij
               dij = -ts1*qi*qj*Rij
               virt(1:3,1) = virt(1:3,1) + dij(1)*Rij
               virt(1:3,2) = virt(1:3,2) + dij(2)*Rij
               virt(1:3,3) = virt(1:3,3) + dij(3)*Rij
            enddo  ! do jsp=1,4
         enddo  ! do isp=1,4
      enddo
   enddo  ! do i=1, Nw
!enddo  ! do j=1, Nw

! print*,'Uelec-1=',Uelec
!print*,'Umon = ' , Umon
!!print*,'Uvdw = ' , Uvdw
!Uelec = 0.5d0 * sum(charge(fO:lH)*phi(fO:lH))
! print*,'Uelec0=',Uelec
!stop

!... self energy
do iat=fO, lH
   phi(iat) = phi(iat) - aewald*TSP*charge(iat) 
enddo
do iw=1, Nwaters
   do isp=1, 2
      iat=iw + (isp-1)*Nwaters
      do jsp=isp+1, 3
         jat=iw + (jsp-1)*Nwaters
         Rij=(R(:,iat)-R(:,jat))
	 Rij=Rij - box*anint( Rij/box )!PBC	
	 dRsq = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
!....
         dd = dsqrt(drsq)
         expon = dexp(-(aewald*aewald*drsq))
         kr0=(erfc(aewald*dd))/dd
         kr1=(kr0 + const_ts1*expon)/drsq
         ch0 = -1.d0/dd
         ch1 = ch0/dRsq
         ts0 = ch0 + kr0
         ts1 = ch1 + kr1
!....
         phi(iat)=phi(iat)+ts0*charge(jat)
         phi(jat)=phi(jat)+ts0*charge(iat)
         Efq(1:3,iat) = Efq(1:3,iat) + ts1*charge(jat)*Rij
         Efq(1:3,jat) = Efq(1:3,jat) - ts1*charge(iat)*Rij
         dij = -ts1*charge(iat)*charge(jat)*Rij
         virt(1:3,1) = virt(1:3,1) + dij(1)*Rij
         virt(1:3,2) = virt(1:3,2) + dij(2)*Rij
         virt(1:3,3) = virt(1:3,3) + dij(3)*Rij
      enddo
   enddo
enddo

!...................................................................
!...............   Reciprocal space  ...............................
!...................................................................
!...
call ewald_pw_std_qq
virt = virt + vir_rec
!...
Uelec = 0.5d0 * sum(charge(fO:lH)*phi(fO:lH))
!print*,'Uelec=',Uelec
do iat=fO, lH
   dR(1:3, iat) = dR(1:3, iat) - charge(iat) * Efq(1:3,iat)
enddo


virt(1,1) = virt(1,1) - Uvdw_lrc
virt(2,2) = virt(2,2) - Uvdw_lrc
virt(3,3) = virt(3,3) - Uvdw_lrc

!..............
dip_mom(1) = sum(charge(fO:lH)*R(1, fO:lH))
dip_mom(2) = sum(charge(fO:lH)*R(2, fO:lH))
dip_mom(3) = sum(charge(fO:lH)*R(3, fO:lH))
do iw=1, Nwaters
   iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters  
   !Rewritten to include PBCs
   rh1O=R(:,iH1)-R(:,iO)
   rh1O = rh1O - box*anint(rh1O/box) !PBC
   
   rh2O=R(:,iH2)-R(:,iO)
   rh2O = rh2O - box*anint(rh2O/box) !PBC

   dip_momI(1:3, iw) = charge(iH1)*rh1O + charge(iH2)*rh2O  

enddo
En = Umon + Uvdw + Uelec + Uvdw_lrc

do iw=1, Nwaters
   iOa=3*iw-2; iH1a=3*iw-1; iH2a=3*iw
   iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters;
   dRR(1:3, (/iOa, ih1a, ih2a/) ) = dR(1:3, (/iO, iH1, iH2/))
enddo

end subroutine pot_spc
