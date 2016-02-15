subroutine ewald_std_qq
use consts
use system_mod
use pot_mod
implicit none
double precision, dimension(3) :: kveck
double precision :: ksq, Akk, p4vAkk, tmp, rQqk, iQqk
integer :: lastkx, lastky
integer :: iat, k, kx, ky,kz
integer :: ndim1, ndim2, ndim3, ntable, nwork, ngrid, npoint, k0, it1, it2, it3, i0, j0, i, j
integer :: nff, nf1, nf2, nf3, k1, k2, k3, m1, m2, m3
double precision :: w, qi, term, t1, t2, t3, pterm, volterm, ar1, ar2, ar3, h1, h2, h3, hsq, erec2
double precision :: expterm, denom, struc2, e, Erec, vterm, dn1, dn2, dn3, de1, de2, de3, s1


call ssincos
lastkx=-1000; lastky=-1000
do k=1, nkvec
   kx = k_hkl(1, k)
   ky = k_hkl(2, k)
   kz = k_hkl(3, k)
   kveck(1:3) = kvec(1:3, k)
   Akk = Ak(k)

!   p4vAkk = p4v*Akk
   p4vAkk = FOURPI*Akk/Volume
   if ( .not. (kx==lastkx .and. ky==lastky)) then
   !... loop over charges 
      do iat=fO, lM
         tmp = sinkx(iat,kx)*cosky(iat,ky) + coskx(iat,kx)*sinky(iat,ky)
         coskxky(iat)= coskx(iat,kx)*cosky(iat,ky) - sinkx(iat,kx)*sinky(iat,ky)
         sinkxky(iat)= tmp
      enddo
   endif
   lastkx=kx;  lastky=ky
   do iat=fO, lM
      tmp = sinkxky(iat)*coskz(iat,kz)+coskxky(iat)*sinkz(iat,kz)
      coskr(iat) = coskxky(iat)*coskz(iat,kz)-sinkxky(iat)*sinkz(iat,kz)
      sinkr(iat) = tmp 
   enddo
   rQqk = sum(charge(fH:lM)*coskr(fH:lM))
   iQqk = sum(charge(fH:lM)*sinkr(fH:lM))
   rQq(k) = rQqk
   iQq(k) = iQqk
   do iat=fO, lM
      tmp = p4vAkk * (coskr(iat)*iQqk - sinkr(iat)*rQqk) 
      Efq(1:3,iat) = Efq(1:3,iat) - tmp*kveck(1:3)
   enddo
   phi(fH:lM) = phi(fH:lM) + p4vAkk * (coskr(fH:lM)*rQqk + sinkr(fH:lM)*iQqk)
enddo

end subroutine ewald_std_qq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ewald_std_dd
use consts
use system_mod
use pot_mod
implicit none
double precision, dimension(3) :: kveck
double precision :: Akk, p4vAkk, tmp, rQqk, iQqk, rQmk, iQmk
integer :: lastkx, lastky
integer :: ii, i, iat, k, kx, ky,kz
integer :: ndim1, ndim2, ndim3, ntable, nwork, ngrid, npoint, k0, it1, it2, it3, i0, j0, j
integer :: nff, nf1, nf2, nf3, k1, k2, k3, m1, m2, m3
double precision :: w, qi, term, t1, t2, t3, pterm, term1,term2,volterm, ar1, ar2, ar3, h1, h2, h3, hsq
double precision :: expterm, denom, struc2, e, Erec, vterm, dn1, dn2, dn3, de1, de2, de3,mu1,mu2,mu3
double precision :: ft1, ft2, ft3, s1, s2, s3
double precision :: du1_dx,du1_dy,du1_dz,  du2_dx,du2_dy,du2_dz,  du3_dx,du3_dy,du3_dz
double precision ::  dphi_du1, dphi_du2 , dphi_du3 

do ii=1, 3
   dipt(fd:ld, ii) = dip(ii, fd:ld)
enddo

lastkx = -1000
lastky = -1000

do k=1, nkvec
   kx = k_hkl(1, k);
   ky = k_hkl(2, k);
   kz = k_hkl(3, k);
   kveck(1:3) = kvec(1:3, k)
   Akk = Ak(k)
!   p4vAkk = p4v*Ak(k)
   p4vAkk = FOURPI*Akk/Volume
   if (.not. (kx==lastkx)) dKkX(fd:ld) = kveck(1)*dipt(fd:ld,1)

   if (.not. (kx==lastkx .and. ky==lastky))  then
      dKkXY(fd:ld) = dKkX(fd:ld) + kveck(2)*dipt(fd:ld,2)
      do i=fd, ld
         tmp=sinkx(i,kx)*cosky(i,ky)+coskx(i,kx)*sinky(i,ky)
         coskxky(i)=coskx(i,kx)*cosky(i,ky)-sinkx(i,kx)*sinky(i,ky)
         sinkxky(i)=tmp
      enddo
   endif

   dKk(fd:ld) = dKkXY(fd:ld) + kveck(3)*dipt(fd:ld,3)
   rQmk = 0.d0
   iQmk = 0.d0
   do i=fd, ld
      tmp=sinkxky(i)*coskz(i,kz) + coskxky(i)*sinkz(i,kz)
      coskr(i)=coskxky(i)*coskz(i,kz) - sinkxky(i)*sinkz(i,kz)
      sinkr(i)=tmp
      rQmk = rQmk - dKk(i)*sinkr(i)
      iQmk = iQmk + dKk(i)*coskr(i)
   enddo
   rQm(k) = rQmk
   iQm(k) = iQmk
   lastkx=kx
   lastky=ky

   do i=fd, ld
      tmp = p4vAkk*(coskr(i)*iQmk-sinkr(i)*rQmk)
      Efd(1:3,i) = Efd(1:3,i) - kveck(1:3) * tmp
   enddo
enddo  ! do while (iblock<=NblocksK)
end subroutine ewald_std_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ewald_std_qd
use consts
use system_mod
use pot_mod
implicit none
double precision, dimension(3) :: kveck, rdipQ, idipQ
double precision, dimension(3,3) :: virtmp, virtmp1
double precision :: ksq, Akk, p4vAkk, tmp, rQqk, iQqk, rQmk, iQmk, rQqmk, iQqmk, coeff2, Ureck
integer :: lastkx, lastky
integer :: ii, i, iat, k, kx, ky,kz,  ix, iy, iz
integer :: ndim1, ndim2, ndim3, ntable, nwork, ngrid, npoint, k0, it1, it2, it3, i0, j0, j
integer :: nff, nf1, nf2, nf3, k1, k2, k3, m1, m2, m3
double precision :: w, qi, term, t1, t2, t3, pterm, term1,term2,volterm, ar1, ar2, ar3, h1, h2, h3, hsq
double precision :: expterm, denom, struc2, e, Erec, vterm, dn1, dn2, dn3, de1, de2, de3,mu1,mu2,mu3
double precision :: ft1, ft2, ft3, ft22, ft23, ft33, s1, s2, s3, s1p
double precision :: du1_dx,du1_dy,du1_dz,  du2_dx,du2_dy,du2_dz,  du3_dx,du3_dy,du3_dz
double precision ::  dphi_du1, dphi_du2 , dphi_du3, de_du1, de_du2, de_du3
double precision :: d2phi_du1du1, d2phi_du1du2, d2phi_du1du3, d2phi_du2du2, d2phi_du2du3, d2phi_du3du3
double precision :: dfx, dfy, dfz, phiiat


sum_rec(fd:ld, 1:3) = 0.d0
do i=1, 3
   dipt(fd:ld, i) = dip(i, fd:ld)
enddo
lastkx=-1000
lastky=-1000
virtmp =0.d0
virtmp1=0.d0
do k=1, nkvec
   kx = k_hkl(1, k);
   ky = k_hkl(2, k);
   kz = k_hkl(3, k);
   kveck(1:3) = kvec(1:3, k)
   Akk = Ak(k)
!   p4vAkk = p4v*Akk
   p4vAkk = FOURPI*Akk/Volume
   if (.not. (kx==lastkx)) dKkX(fd:ld) = kveck(1)*dipt(fd:ld,1)

   if (.not. (kx==lastkx .and. ky==lastky))  then
      dKkXY(fd:ld) = dKkX(fd:ld) + kveck(2)*dipt(fd:ld,2)
      do i=fO, lM
         tmp=sinkx(i,kx)*cosky(i,ky)+coskx(i,kx)*sinky(i,ky)
         coskxky(i)=coskx(i,kx)*cosky(i,ky)-sinkx(i,kx)*sinky(i,ky)
         sinkxky(i)=tmp
      enddo
   endif

   dKk(fd:ld) = dKkXY(fd:ld) + kveck(3)*dipt(fd:ld,3)
   do i=fO, lM
      tmp=sinkxky(i)*coskz(i,kz) + coskxky(i)*sinkz(i,kz)
      coskr(i)=coskxky(i)*coskz(i,kz) - sinkxky(i)*sinkz(i,kz)
      sinkr(i)=tmp
   enddo
   lastkx=kx
   lastky=ky

   rQmk = rQm(k);           iQmk  = iQm(k)
   rQqk = rQq(k);           iQqk  = iQq(k)
   rQqmk = rQq(k) + rQm(k); iQqmk = iQq(k) + iQm(k)

   do iat=fH, lM    ! charges
      tmp = charge(iat)*p4vAkk*(coskr(iat)*iQmk-sinkr(iat)*rQmk)
      phi(iat) = phi(iat) + p4vAkk*(coskr(iat)*rQmk + sinkr(iat)*iQmk)
      dR(1:3,iat) = dR(1:3,iat) + tmp*kveck(1:3)
   enddo

   do iat=fd, ld    ! dipoles
      tmp =  dkk(iat)*p4vAkk*(coskr(iat)*rQqmk+sinkr(iat)*iQqmk)
      dR(1:3,iat) = dR(1:3,iat) - tmp*kveck(1:3)
   enddo
   !... virial
   ksq=kveck(1)*kveck(1) + kveck(2)*kveck(2) + kveck(3)*kveck(3)
   coeff2 = 2.d0 / ksq + 0.5d0  / (aewald*aewald)
   Ureck = Akk*(rQqmk**2 + iQqmk**2)
   virtmp(1,1) = virtmp(1,1) + Ureck - Ureck * coeff2*kveck(1)*kveck(1)
   virtmp(2,1) = virtmp(2,1)         - Ureck * coeff2*kveck(1)*kveck(2)
   virtmp(3,1) = virtmp(3,1)         - Ureck * coeff2*kveck(1)*kveck(3)
   virtmp(2,2) = virtmp(2,2) + Ureck - Ureck * coeff2*kveck(2)*kveck(2)
   virtmp(3,2) = virtmp(3,2)         - Ureck * coeff2*kveck(2)*kveck(3)
   virtmp(3,3) = virtmp(3,3) + Ureck - Ureck * coeff2*kveck(3)*kveck(3)
   virtmp(1,2) = virtmp(1,2)         - Ureck * coeff2*kveck(1)*kveck(2)
   virtmp(1,3) = virtmp(1,3)         - Ureck * coeff2*kveck(1)*kveck(3)
   virtmp(2,3) = virtmp(2,3)         - Ureck * coeff2*kveck(2)*kveck(3)

   rdipQ(1) = -sum(sinkr(fd:ld)*dipt(fd:ld, 1) )
   rdipQ(2) = -sum(sinkr(fd:ld)*dipt(fd:ld, 2) )
   rdipQ(3) = -sum(sinkr(fd:ld)*dipt(fd:ld, 3) )
   idipQ(1) =  sum(coskr(fd:ld)*dipt(fd:ld, 1) )
   idipQ(2) =  sum(coskr(fd:ld)*dipt(fd:ld, 2) )
   idipQ(3) =  sum(coskr(fd:ld)*dipt(fd:ld, 3) )

   do ix=1,3
      tmp = 2.d0*Akk*kveck(ix)
      do iy=1,3
         virtmp1(ix,iy)=virtmp1(ix,iy) + tmp*(rQqmk*rdipQ(iy)+iQqmk*idipQ(iy))
      enddo
   enddo
enddo
virtmp  = -0.5d0*virtmp*FOURPI/Volume
virtmp1 = -0.5d0*virtmp1*FOURPI/Volume
vir_rec = vir_rec + virtmp + virtmp1
   
end subroutine ewald_std_qd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ewald_pw_std_qq   ! for pairwise additive
use consts
use system_mod
use pot_mod
implicit none
double precision, dimension(3) :: kveck
double precision :: ksq, Akk, p4vAkk, tmp, rQqk, iQqk
integer :: lastkx, lastky
integer :: iat, k, kx, ky,kz
integer :: ndim1, ndim2, ndim3, ntable, nwork, ngrid, npoint, k0, it1, it2, it3, i0, j0, i, j
integer :: nff, nf1, nf2, nf3, k1, k2, k3, m1, m2, m3
double precision :: w, qi, term, t1, t2, t3, pterm, volterm, ar1, ar2, ar3, h1, h2, h3, hsq, erec2
double precision :: expterm, denom, struc2, e, Erec, vterm, dn1, dn2, dn3, de1, de2, de3, s1, Ureck, coeff2

vir_rec = 0.d0
call ssincos
lastkx=-1000; lastky=-1000
do k=1, nkvec
   kx = k_hkl(1, k)
   ky = k_hkl(2, k)
   kz = k_hkl(3, k)
   kveck(1:3) = kvec(1:3, k)
   Akk = Ak(k)
!   p4vAkk = p4v*Akk
   p4vAkk = FOURPI*Akk/Volume
   if ( .not. (kx==lastkx .and. ky==lastky)) then
   !... loop over charges 
      do iat=fO, lH
         tmp = sinkx(iat,kx)*cosky(iat,ky) + coskx(iat,kx)*sinky(iat,ky)
         coskxky(iat)= coskx(iat,kx)*cosky(iat,ky) - sinkx(iat,kx)*sinky(iat,ky)
         sinkxky(iat)= tmp
      enddo
   endif
   lastkx=kx;  lastky=ky
   do iat=fO, lH
      tmp = sinkxky(iat)*coskz(iat,kz)+coskxky(iat)*sinkz(iat,kz)
      coskr(iat) = coskxky(iat)*coskz(iat,kz)-sinkxky(iat)*sinkz(iat,kz)
      sinkr(iat) = tmp 
   enddo
   rQqk = sum(charge(fO:lH)*coskr(fO:lH))
   iQqk = sum(charge(fO:lH)*sinkr(fO:lH))
   do iat=fO, lH
      tmp = p4vAkk * (coskr(iat)*iQqk - sinkr(iat)*rQqk) 
      Efq(1:3,iat) = Efq(1:3,iat) - tmp*kveck(1:3)
   enddo
   phi(fO:lH) = phi(fO:lH) + p4vAkk * (coskr(fO:lH)*rQqk + sinkr(fO:lH)*iQqk)
   !... virial
   ksq=kveck(1)*kveck(1) + kveck(2)*kveck(2) + kveck(3)*kveck(3)
   coeff2 = 2.d0 / ksq + 0.5d0  / (aewald*aewald)
    Ureck= -0.5d0*p4VAkk*(rQqk**2 + iQqk**2)
   vir_rec(1,1) = vir_rec(1,1) + Ureck - Ureck * coeff2*kveck(1)*kveck(1)
   vir_rec(2,1) = vir_rec(2,1)         - Ureck * coeff2*kveck(1)*kveck(2)
   vir_rec(3,1) = vir_rec(3,1)         - Ureck * coeff2*kveck(1)*kveck(3)
   vir_rec(2,2) = vir_rec(2,2) + Ureck - Ureck * coeff2*kveck(2)*kveck(2)
   vir_rec(3,2) = vir_rec(3,2)         - Ureck * coeff2*kveck(2)*kveck(3)
   vir_rec(3,3) = vir_rec(3,3) + Ureck - Ureck * coeff2*kveck(3)*kveck(3)
   vir_rec(1,2) = vir_rec(1,2)         - Ureck * coeff2*kveck(1)*kveck(2)
   vir_rec(1,3) = vir_rec(1,3)         - Ureck * coeff2*kveck(1)*kveck(3)
   vir_rec(2,3) = vir_rec(2,3)         - Ureck * coeff2*kveck(2)*kveck(3)
enddo

end subroutine ewald_pw_std_qq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ewald_set(first)
use pot_mod
use system_mod
use consts
implicit none
integer :: grp
integer :: nmaxkvec, kx, ky, kz, ky1, kz1
double precision, dimension(3) :: kv
double precision :: PP, kvec_cut, kvec_cutSQ, factor, kvsq, kvxy
double precision :: fac_ew2
integer :: i, k, nf1
integer :: ifft1, ifft2, ifft3, ntable, itmp, my_3S, my_3E, my_2S, my_2E
integer, parameter :: maxpower=50
logical :: FIRST

if (FIRST) then 
!/////////////////////////////////////////////////////////!
!////////////  EWALD SET-UP///////////////////////////////!
!/////////////////////////////////////////////////////////!
PP = -dlog(eps_ewald)
aewald = dsqrt(PP)/Rc

   kvec_cut = 2.d0*aewald*dsqrt(PP)
   kvec_cutSQ = kvec_cut**2
   hmax = floor(kvec_cut/TWOPI*box(1))
   kmax = floor(kvec_cut/TWOPI*box(2))
   lmax = floor(kvec_cut/TWOPI*box(3))
   fac_ew2 = -0.25d0/(aewald*aewald)

      !print*,"Ewald parameter : ", aewald
     ! print*,'hmax=',hmax

   nmaxkvec = 4*(hmax+1)*(kmax+1)*(lmax+1)

   allocate(kvec(3, nmaxkvec))
   allocate(Ak(nmaxkvec))
   allocate(k_hkl(3, nmaxkvec))
   allocate(coskx(1:NatomsM, 0:hmax))
   allocate(cosky(1:NatomsM, -kmax:kmax))
   allocate(coskz(1:NatomsM, -lmax:lmax))
   allocate(sinkx(1:NatomsM, 0:hmax))
   allocate(sinky(1:NatomsM, -kmax:kmax))
   allocate(sinkz(1:NatomsM, -lmax:lmax))
   allocate(coskxky(1:NatomsM))
   allocate(sinkxky(1:NatomsM))
   allocate(coskr(1:NatomsM))
   allocate(sinkr(1:NatomsM))
   allocate(rQm(nmaxkvec))
   allocate(iQm(nmaxkvec))
   allocate(rQq(nmaxkvec))
   allocate(iQq(nmaxkvec))

   nkvec = 0
   factor=2.d0
endif 

   do kx=0, hmax
      kv(1) = TWOPI*boxi(1)*dble(kx)
      if (kx==0) then
         ky1 = 0
      else
         ky1 = -kmax
      endif

      do ky=ky1, kmax
         kv(2) = TWOPI*boxi(2)*dble(ky)
          if (kx==0 .and. ky==0) then
             kz1=1
          else
             kz1=-lmax
          endif

         kvxy = kv(1)**2 + kv(2)**2
         do kz=kz1, lmax
            kv(3) = TWOPI*boxi(3)*dble(kz)
            kvsq = kvxy + kv(3)**2
            if (kvsq<kvec_cutSQ .and. kvsq> 1.d-16) then
               nkvec = nkvec+1
               k_hkl(1:3, nkvec) = (/kx, ky, kz/)
               kvec(1:3, nkvec) = kv(1:3)
               Ak(nkvec) = factor*dexp(fac_ew2*kvsq)/kvsq
            endif
         enddo
      enddo
   enddo
end subroutine ewald_set
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine ewald_scl_kvec(grp)
use pot_mod
use system_mod
use consts
implicit none
integer :: grp
integer :: ik, kx, ky, kz, ky1, kz1
double precision :: factor, kvsq, fac_ew2, kvxy, kvec_cut , kvec_cutSQ, PP
double precision, dimension(3) :: kv


fac_ew2 = -0.25d0/(aewald*aewald)
factor = 2.d0
do ik=1, nkvec
   kvec(1:3,ik) = TWOPI*dble(k_hkl(1:3,ik))*boxi(1:3)
   kvsq = kvec(1,ik)*kvec(1,ik) + kvec(2,ik)*kvec(2,ik) + kvec(3,ik)*kvec(3,ik)
   Ak(ik) = factor*dexp(fac_ew2*kvsq)/kvsq
enddo


end subroutine ewald_scl_kvec
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine ssincos
use consts
use system_mod
use pot_mod
implicit none
integer :: i, kx, ky, kz
double precision :: vx, vy, vz
double precision, dimension(3) :: Ri

coskx(:, 0) = 1.d0; cosky(:, 0)=1.d0; coskz(:, 0)=1.d0
sinkx(:, 0) = 0.d0; sinky(:, 0)=0.d0; sinkz(:, 0)=0.d0
coskx(:, 1)=dcos(TWOPI*boxi(1)*R(1,:));
sinkx(:, 1)=dsin(TWOPI*boxi(1)*R(1,:));
cosky(:, 1)=dcos(TWOPI*boxi(2)*R(2,:));cosky(:,-1)=cosky(:, 1)
sinky(:, 1)=dsin(TWOPI*boxi(2)*R(2,:));sinky(:,-1)=-sinky(:,1)
coskz(:, 1)=dcos(TWOPI*boxi(3)*R(3,:));coskz(:,-1)=coskz(:, 1)
sinkz(:, 1)=dsin(TWOPI*boxi(3)*R(3,:));sinkz(:,-1)=-sinkz(:,1)
DO  KX = 2, hmax
   coskx(:, KX) = coskx(:, kx-1)*coskx(:,1) - sinkx(:, kx-1)*sinkx(:,1)
   sinkx(:, KX) = coskx(:, kx-1)*sinkx(:,1) + sinkx(:, kx-1)*coskx(:,1)
ENDDO
DO  ky = 2, kmax
   cosky(:, ky) = cosky(:, ky-1)*cosky(:,1) - sinky(:, ky-1)*sinky(:,1)
   sinky(:, ky) = cosky(:, ky-1)*sinky(:,1) + sinky(:, ky-1)*cosky(:,1)
   cosky(:,-ky) = cosky(:, ky)
   sinky(:,-ky) =-sinky(:, ky)
ENDDO
DO  kz = 2, lmax
   coskz(:, kz) = coskz(:, kz-1)*coskz(:,1) - sinkz(:, kz-1)*sinkz(:,1)
   sinkz(:, kz) = coskz(:, kz-1)*sinkz(:,1) + sinkz(:, kz-1)*coskz(:,1)
   coskz(:,-kz) = coskz(:, kz)
   sinkz(:,-kz) =-sinkz(:, kz)
ENDDO

end subroutine ssincos
