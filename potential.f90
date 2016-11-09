!---------------------------------------------------------------------!
!-----------------Call correct potential -----------------------------
!---------------------------------------------------------------------!
subroutine potential(RR, RRc, Upot, dRR, virt, virialc, dip_momI, Edip_mom, chg, t, BAROSTAT)
use consts
use system_mod
use pot_mod
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR, RRc
double precision, intent(out) :: Upot, virialc
double precision, dimension(3, Natoms), intent(out) :: dRR
double precision, dimension(3, 3), intent(out) :: virt
double precision, dimension(Natoms), intent(out)  :: chg
double precision, dimension(3, NWaters), intent(out)  ::  dip_momI, Edip_mom
double precision, dimension(3) :: dip_mom
double precision, dimension(3) :: tmp_der
integer :: i
integer, intent(in) :: t
logical, intent(in) :: BAROSTAT

!All the stuff that depends on volume needs to be rescaled. 
 	p4V = FOURPI/volume

!scale VdW long range correction due to box size change and add correction
!multiplying by a correction factor is slightly more efficient than recalculating the entire Uvdw_lrc term each timestep
	Uvdw_lrc = Uvdw_lrc0*(volume_init/volume)

!If the volume is changing than the Ewald k-vectors have to be reset every timestep
!if (BAROSTAT) call ewald_set(.false.)

if (pot_model==2 .or. pot_model==3) then
    call pot_ttm(RR, RRc, Upot, dRR, virt, virialc, dip_momI, Edip_mom, chg,t)

else if (pot_model==4 .or. pot_model==5) then
    call pot_spc(RR, Upot, dRR, virt, dip_momI, chg)
endif

end subroutine potential

!---------------------------------------------------------------------!
!-----------------Initialize potential variables ---------------------
!---------------------------------------------------------------------!
subroutine init_pot
use consts
use system_mod
use nasa_mod
use pot_mod
implicit none
double precision :: polfacO, polfacH, polfacM
integer :: iw, io, ih1, ih2


if (pot_model==2 .or. pot_model==3) then
   polar_model = .true.
else
   polar_model = .false.
endif

if (allocated(R)) then
   deallocate(R, dR, charge, Efq, phi)
   if (polar_model) deallocate(Efd, dip, olddip, grdq)
endif

if (polar_model) then
   Nwaters = Natoms/3
   NatomsM = 4*Nwaters
   fO=1;   lO=Nwaters
   fH=Nwaters+1; lH=3*Nwaters
   fM=3*Nwaters+1; lM=4*Nwaters
   fO3=3*fO-2; lO3=3*lO
   fH3=3*fH-2; lH3=3*lH
   fM3=3*fM-2; lM3=3*lM
else
   Nwaters = Natoms/3
   NatomsM = 3*Nwaters
   fO=1;   lO=Nwaters
   fH=Nwaters+1; lH=3*Nwaters
   fM=Nwaters+1; lM=3*Nwaters
   fO3=3*fO-2; lO3=3*lO
   fH3=3*fH-2; lH3=3*lH
   fM3=3*fH-2; lM3=3*lH
endif

if (pot_model==2) then
  ! print*,'>>>>>>>>> THE MODEL IS TTM2.1-F <<<<<<<<<<<<<<'
   fd=fO; fd3=fO3; ld=lH; ld3=lH3
   vdwA = ttm21f_vdwA; vdwB = ttm21f_vdwB; vdwC = ttm21f_vdwC; 
   vdwD = ttm21f_vdwD; vdwE = ttm21f_vdwE; 
   polfacO=ttm21f_polfacO; polfacH=ttm21f_polfacH; polfacM=ttm21f_polfacM;
   polarO=ttm21f_polarO; polarH=ttm21f_polarH; polarM=ttm21f_polarM;
   gammaM=ttm21f_gammaM; aDD = ttm21f_aDD; aCCaCD=ttm21f_aCCaCD
   dms_param1 = ttm21f_dms_param1; dms_param2 = ttm21f_dms_param2;
   dms_param3 = ttm21f_dms_param3
   fdI = 1; ldI=3
else if (pot_model==3) then
  ! print*,'>>>>>>>>> THE MODEL IS TTM3-F <<<<<<<<<<<<<<'
   fd=fM; fd3=fM3; ld=lM; ld3=lM3
   vdwA=ttm3f_vdwA; vdwB=ttm3f_vdwB; vdwC=ttm3f_vdwC; vdwD=ttm3f_vdwD; vdwE=ttm3f_vdwE;
   polfacO=ttm3f_polfacO; polfacH=ttm3f_polfacH; polfacM=ttm3f_polfacM;
   polarO=ttm3f_polarO; polarH=ttm3f_polarH; polarM=ttm3f_polarM;
   gammaM=ttm3f_gammaM; aDD = ttm3f_aDD; aCCaCD=ttm3f_aCCaCD;
   dms_param1 = ttm3f_dms_param1; dms_param2 = ttm3f_dms_param2;
   dms_param3 = ttm3f_dms_param3
   fdI = 4; ldI=4
        
   !GROMACS VdW shift parameters. see GROMACS manual Ch. 4
   if (Rc .eq. rc1) then
	shiftA = 0
	shiftB = 0 
   else
	shiftA = - (10*Rc - 7*rc1 )/(Rc**8 * (Rc - rc1)**2) 
   	shiftB =   (9*Rc - 7*rc1 )/ (Rc**8 * (Rc - rc1)**3)
   endif


   !Coloumb shift parameters (not used & not recommended!)
   !shiftAc = -(5*Rc - 2*rc1 )/(Rc**3 * (Rc - rc1)**2 )
   !shiftBc =  (4*Rc - 2*rc1 )/(Rc**3 * (Rc - rc1)**3 )


else if (pot_model==4) then
   !print*,'>>>>>>>>> THE MODEL IS qSPC/Fw (Voth) <<<<<<<<<<<<<<'
   vdwA=qspcfw_vdwA; vdwB=qspcfw_vdwB; vdwC=qspcfw_vdwC; vdwD=qspcfw_vdwD; vdwE=qspcfw_vdwE;
   qO = qspcfw_qO; qH=qspcfw_qH
   roh0 = qspcfw_roh0; th0 = qspcfw_th0
   harmkb = qspcfw_harmkb; harmka = qspcfw_harmka

   !vdwA=629326.57249877361701443513d0
   !vdwB=0.d0
   !vdwC=-625.50206521141523463400d0
   !vdwD=0.d0; 
   !vdwE=0.d0;
   !qO = -0.84d0; 
   !qH=0.42d0
   !roh0 = 1.012d0; 
   !th0 = 1.95476876223364912614d0 ! =112.0 deg
   !harmkb = 529.581d0;
   !harmka = 37.95d0

else if (pot_model==5) then
  ! print*,'>>>>>>>>> THE MODEL IS  SPC/F JCP_106_2400 > >>>>>>>>>>>>> <<<<<<<<<<<<<<'
   vdwA=spcf_vdwA; vdwB=spcf_vdwB; vdwC=spcf_vdwC; vdwD=spcf_vdwD; vdwE=spcf_vdwE;
   qO = spcf_qO; qH=spcf_qH
   roh0 = spcf_roh0; th0 = spcf_th0
endif

Rc_sq = Rc*Rc
volume = box(1)*box(2)*box(3)
boxi = 1.d0 / box
p4V = FOURPI/volume !keeping this variable updated is problematic when the volume is changing.  

Uvdw_lrc= TWOPI*dble(Nwaters*Nwaters)/volume * &
            (vdwA/(9.d0*(Rc**9)) + vdwB/(7.d0*(Rc**7)) +vdwC/(3.d0*(Rc**3)))
Uvdw_lrc0 = Uvdw_lrc

!write(*,*) Uvdw_lrc, "= long range correction"
if (polar_model) then
   polfac(1:4)=(/ polfacO**(1.d0/6.d0), polfacH**(1.d0/6.d0), polfacH**(1.d0/6.d0), polfacM**(1.d0/6.d0) /)
   polfac(1:4)=1.d0/polfac(1:4)**3
endif



allocate(neigh_list(2, Nwaters*Nwaters/2))
allocate(R(3, NatomsM))
allocate(dR(3, NatomsM))
allocate(charge(NatomsM))
allocate(Efq(3,NatomsM))
allocate(phi(NatomsM))
allocate(dKk(NatomsM))
allocate(dKkX(NatomsM))
allocate(dKkXY(NatomsM))
allocate(sum_rec(fd:ld, 3))


if (polar_model) then
   allocate(Efd(3,NatomsM))
   allocate(dip(3,fd:ld))
   allocate(dipt(fd:ld, 3))
   allocate(olddip(3,fd:ld))
   allocate(tx_dip(3,fd:ld, 4))  
   allocate(grdq(Nwaters,3,3,3))
   tx_dip = 0 
endif

if (pot_model==4 .or. pot_model==5) then
   do iw=1, Nwaters
!      iO=3*iw-2; ih1=3*iw-1; ih2=3*iw
      iO=fO+iw-1; iH1=iO+Nwaters; iH2=iO+2*Nwaters
      charge( (/iO, ih1, ih2/) ) = (/qO, qh, qh/) * CHARGECON
   enddo
endif


call ewald_set(.true.) !initial Ewald set up 
const_ts1 = 2.d0*aewald/SQPI
const_ts2 = 4.d0*aewald*aewald*aewald/(3.d0*SQPI)
const_ts3 = 8.d0*(aewald**5)/(15.d0*SQPI)
end subroutine init_pot
