!----------------------------------------------------
! Dipole moments from TTM3F
! scaled down version of pot_ttm.f90 that only calculates dipole moments
!----------------------------------------------------
subroutine dip_ttm(RR, dip_momI, Edip_mom, t)
use consts
use system_mod
use pot_mod
use neigh_mod
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR 
 double precision, dimension(3, NWaters), intent(out) :: dip_momI, Edip_mom
double precision, dimension(Natoms), intent(out) :: chg
integer :: itmp, i, is, j, js, iw, jw, iO, ih1, ih2, iOa, iH1a, iH2a, iM, jO, jH1, jH2, jM
integer :: iat, jat, isp, jsp, jsp0, i3, j3, ii, jj, kk, ix, iy
double precision :: ts0, ts1, ts2, ts3, ts1DD, ts2DD, ts3DD, pol, qi, qj, Ureck, e1
double precision :: tmp, R2i, R4i, R6i, R8i, uij, dRsq, dRij, qdqd_sq, ksq, polfacI
double precision :: dri, drsqi, expon, er, dd, ch1, ch2, ch3, sr1, sr2, sr3, exp1, a_pol12, drijsq
double precision :: Uind_init
double precision, dimension(3) :: qdqd,q3, Ri,Rij, dij,di,dj, roh1,roh2,rhh, Rij0,rh1m,rh2m,rmO, Rcij
double precision, dimension(3, 3) :: r1, dr1, vij, dd3, arr33, dgrad
double precision, dimension(3,3,3) :: dq3 
double precision, dimension(4) :: dpcf
! double precision :: factor, deltadip, stath 
integer :: iter
integer, intent(in) :: t
logical :: EWALD = .true.
type(t_neigh) :: neigh 

dip = 0.d0
olddip = 0.d0


!converts OHHOHH... to OOOO...HHHH...HHHH...MMMM
do iw=1, Nwaters
   iOa=3*iw-2; iH1a=3*iw-1; iH2a=3*iw
   iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters; iM=iO+3*Nwaters
   R(1:3, (/iO, iH1, iH2/)) = RR(1:3, (/iOa, ih1a, ih2a/) )! * boxi
   RRRc(1:3, (/iO, iH1, iH2/)) = RRc(1:3, (/iOa, ih1a, ih2a/) )! * boxi
  
!.... determine M-site positions for R
   roh1 = R(1:3, iH1) - R(1:3,iO)
   roh1 = roh1 - box*anint(roh1*boxi)!PBC

   roh2 = R(1:3, iH2) - R(1:3,iO)
   roh2 = roh2 - box*anint(roh2*boxi)!PBC

   R(1:3, iM) = 0.5d0*gammaM*( roh1(:) + roh2(:) ) + R(:,iO)
enddo


tmp = 0.5d0*gammaM/(1.d0-gammaM)
charge = 0.d0
grdq = 0.d0
dR = 0.d0
do iw=1, Nwaters
   iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters; iM=iO+3*Nwaters
   r1(1:3, 1:3) = R(1:3, (/iO, ih1, ih2/) )


   call dms_nasa(r1, q3, dq3,box,boxi)

   q3 = q3*CHARGECON
   dq3 = dq3*CHARGECON

   chg(3*iw-2:3*iw) = q3

   charge(iO) = 0.d0
   charge(iH1) = q3(2) + tmp*(q3(2)+q3(3))
   charge(iH2) = q3(3) + tmp*(q3(2)+q3(3))
   charge(iM ) = q3(1) / (1.d0-gammaM)
 
enddo
 



!........................................................................................!
!.......... Initial guess of induced dipoles ............................................!
!........................................................................................!
predict_step = t
!write(*,*) "in force calc, t = ", t
if ( (predict_step .gt. 4) .and. guess_initdip) then
	!predictor 
	call calc_dpcf(predict_step, dpcf)
	! 	write(*,*) dpcf
	   dip(1:3,fd:ld) = dpcf(1)*tx_dip(1:3,fd:ld, 1) + dpcf(2)*tx_dip(1:3,fd:ld, 2) &
                   + dpcf(3)*tx_dip(1:3,fd:ld, 3) + dpcf(4)*tx_dip(1:3,fd:ld, 4)
else
	!initial guess (t = 1)
	if (pot_model==2) then
		dip(1:3,fO:lO) = polarO * Efq(1:3,fO:lO)
 	        dip(1:3,fH:lH) = polarH * Efq(1:3,fH:lH)
	else
		dip(1:3,fM:lM) = polarM * Efq(1:3,fM:lM)
	endif
endif
Uind_init = 0.0d0
do i=fd, ld
   Uind_init = Uind_init + dip(1,i)*Efq(1,i) + dip(2,i)*Efq(2,i) + dip(3,i)*Efq(3,i)
enddo
Uind_init = -0.5d0*Uind_init
!........................................................................................!
!.......... Calculation ind-dipoles via iterations ......................................!
!........................................................................................!
stath = DEBYE/CHARGECON/dsqrt(dble(Natoms))
factor = 4.d0*(aewald**3)/(3.d0*dsqrt(PI))
do iter=1, polar_maxiter
   Efd = 0.d0

   do i=1, Nwaters
      iO = i
      call find_neigh(iO, R, Neigh)
      do j=1, neigh % N
         jO = neigh % j(j)
         Rij = neigh % Rij(1:3, j)
         dRsq = neigh % R2(j)
         Rij0 = Rij - R(1:3, iO) + R(1:3, jO)!????????????????????
	 Rij0 = Rij0 - box*anint(Rij0*boxi) !PBC
         do isp=fdI, ldI
            iat = iO + (isp-1)*Nwaters
            polfacI = aDD*polfac(isp)
            di = dip(1:3, iat)
            Ri = R(1:3, iat) + Rij0
            do jsp=fdI, ldI
               jat = jO + (jsp-1)*Nwaters
               dj = dip(1:3, jat)
               Rij = Ri - R(1:3, jat)
		Rij = Rij - box*anint(Rij*boxi) !PBC
               dRsq = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
               call PBCsmear2(dRsq, polfacI*polfac(jsp), aewald, ts1DD, ts2DD)
               call get_dd3(Rij, ts1DD, ts2DD, dd3)
               Efd(1,iat) = Efd(1,iat)+dd3(1,1)*dj(1)+dd3(2,1)*dj(2)+dd3(3,1)*dj(3)
               Efd(2,iat) = Efd(2,iat)+dd3(1,2)*dj(1)+dd3(2,2)*dj(2)+dd3(3,2)*dj(3)
               Efd(3,iat) = Efd(3,iat)+dd3(1,3)*dj(1)+dd3(2,3)*dj(2)+dd3(3,3)*dj(3)
               Efd(1,jat) = Efd(1,jat)+dd3(1,1)*di(1)+dd3(2,1)*di(2)+dd3(3,1)*di(3)
               Efd(2,jat) = Efd(2,jat)+dd3(1,2)*di(1)+dd3(2,2)*di(2)+dd3(3,2)*di(3)
               Efd(3,jat) = Efd(3,jat)+dd3(1,3)*di(1)+dd3(2,3)*di(2)+dd3(3,3)*di(3)
            enddo
         enddo
      enddo
   enddo
   !... self interactions
   do iw=1, Nwaters
      do isp=fdI, ldI-1
         iat=iw + (isp-1)*Nwaters
         di = dip(1:3, iat)
         polfacI = aDD*polfac(isp)
         do jsp=isp+1, ldI
            jat=iw + (jsp-1)*Nwaters
            dj = dip(1:3, jat)
            Rij=(R(:,iat)-R(:,jat))
	    Rij = Rij - box*anint(Rij*boxi) !PBC
	    dRsq=Rij(1)*Rij(1)+Rij(2)*Rij(2) + Rij(3)*Rij(3)
            call PBCsmear2(dRsq, polfacI*polfac(jsp), aewald, ts1DD, ts2DD)
            call get_dd3(Rij, ts1DD, ts2DD, dd3)
            Efd(1,iat) = Efd(1,iat)+dd3(1,1)*dj(1)+dd3(2,1)*dj(2)+dd3(3,1)*dj(3)
            Efd(2,iat) = Efd(2,iat)+dd3(1,2)*dj(1)+dd3(2,2)*dj(2)+dd3(3,2)*dj(3)
            Efd(3,iat) = Efd(3,iat)+dd3(1,3)*dj(1)+dd3(2,3)*dj(2)+dd3(3,3)*dj(3)
            Efd(1,jat) = Efd(1,jat)+dd3(1,1)*di(1)+dd3(2,1)*di(2)+dd3(3,1)*di(3)
            Efd(2,jat) = Efd(2,jat)+dd3(1,2)*di(1)+dd3(2,2)*di(2)+dd3(3,2)*di(3)
            Efd(3,jat) = Efd(3,jat)+dd3(1,3)*di(1)+dd3(2,3)*di(2)+dd3(3,3)*di(3)
         enddo
      enddo
   enddo
   ! ................. reciprocal space
  if (EWALD) then 
   call ewald_std_dd
  endif 
   !.............
   Efd(1:3, fd:ld) = Efd(1:3, fd:ld) + Efq(1:3,fd:ld) + factor*dip(1:3,fd:ld)
   !.............. calculate a new induced dipole
   if (pot_model==2) then
	dip(1:3,fO:lO) = (1.d0-polar_sor)*dip(1:3,fO:lO) + polar_sor*polarO*Efd(1:3,fO:lO)
        dip(1:3,fH:lH) = (1.d0-polar_sor)*dip(1:3,fH:lH) + polar_sor*polarH*Efd(1:3,fH:lH)
   else
	dip(1:3,fM:lM) = (1.d0-polar_sor)*dip(1:3,fM:lM) + polar_sor*polarM*Efd(1:3,fM:lM)
   endif
   !..... check for convergence.............
   Uind = -0.5d0*sum(dip(1:3,fd:ld)*Efq(1:3,fd:ld));
   deltadip = stath*dsqrt(sum( (dip-olddip)**2 ))
   if (deltadip<polar_eps) then
goto 100
   endif
olddip = dip
enddo


100 Continue
Uind = -0.5d0*sum(dip(1:3,fd:ld)*Efq(1:3,fd:ld));
itmp = mod(predict_step, 4)+1
if (guess_initdip) tx_dip(1:3, fd:ld, itmp) = dip(1:3, fd:ld)


if (print_dipiters) write(*,'(a,f14.8,2x,a,f14.8,2x,a,i3)') &
    &"Uind(init)=",Uind_init,"Uind(fin)=",Uind,"#of iterations = ",iter

    
   
!calculate final dipole moments
do iw=1, Nwaters
   iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters; iM=iO+3*Nwaters
   !Rewritten to include PBCs
   rh1m=R(:,iH1)-R(:,iM)
   rh1m = rh1m - box*anint(rh1m*boxi) !PBC

   rh2m=R(:,iH2)-R(:,iM)
   rh2m = rh2m - box*anint(rh2m*boxi) !PBC

   rmO =R(:,iM) - R(:,iO)
   rmO = rmO - box*anint(rmO*boxi) !PBC
   
   dip_momI(1:3, iw) = charge(iH1)*rh1m + charge(iH2)*rh2m
 
   !dip_momI(1:3, iw) = charge(iH1)*rh1m + charge(iH2)*rh2m + ( charge(iH1) + charge(iM) + charge(iH2) )*R(:,iM)
   !dip_momI(1:3, iw) = charge(iH1)*R(1:3, iH1) + charge(iH2)*R(1:3, iH2) + charge(iM)*R(1:3, iM)

   do isp=fdI, ldI
       iat= iw + (isp-1)*Nwaters
       dip_momI(1:3, iw) = dip(1:3, iat) + dip_momI(1:3, iw)
       Edip_mom(1:3, iw) = dip(1:3, iat) !sum of polarization dipole(s) for each atom in molecule
   enddo
enddo
end subroutine dip_ttm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_dd3(R, ts1, ts2, dd3)
implicit none
double precision, dimension(3) :: R
double precision :: ts1, ts2
double precision, dimension(3,3) :: dd3

dd3(1,1) = 3.d0*ts2*R(1)*R(1) - ts1
dd3(2,2) = 3.d0*ts2*R(2)*R(2) - ts1
dd3(3,3) = 3.d0*ts2*R(3)*R(3) - ts1
dd3(1,2) = 3.d0*ts2*R(1)*R(2)
dd3(1,3) = 3.d0*ts2*R(1)*R(3)
dd3(2,3) = 3.d0*ts2*R(2)*R(3)
dd3(2,1) = dd3(1,2); dd3(3,1) = dd3(1,3); dd3(3,2) = dd3(2,3)

end subroutine get_dd3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_dpcf(md_step, dpcf)
implicit none
integer, intent(in) :: md_step
double precision, dimension(4), intent(out) :: dpcf
integer :: itmp

if (md_step<4) then
if (md_step==1) then
dpcf(1:4) = (/1.d0, 0.d0, 0.d0, 0.d0/)
   else if (md_step==2) then
dpcf(1:4) = (/-1.d0, 2.d0, 0.d0, 0.d0/)
   else if (md_step==3) then
dpcf(1:4) = (/1.d0, -3.d0, 3.d0, 0.d0/)
   endif
else
itmp = mod(md_step, 4)
   if (itmp==0) then
dpcf(1:4) = (/-1.d0, 4.d0, -6.d0, 4.d0/)
   else if (itmp==1) then
dpcf(1:4) = (/ 4.d0, -1.d0, 4.d0, -6.d0/)
   else if (itmp==2) then
dpcf(1:4) = (/-6.d0, 4.d0, -1.d0, 4.d0 /)
   else if (itmp==3) then
dpcf(1:4) = (/ 4.d0, -6.d0, 4.d0, -1.d0 /)
   endif
endif
end subroutine calc_dpcf
