module pot_mod
implicit none

!----------------------------------------------------------------------------!
! Parameters of the TTM21F potential                                         !
!----------------------------------------------------------------------------!
double precision, parameter :: ttm21f_vdwA=-1.329565985D+6,ttm21f_vdwB=3.632560798D+5,  &
                               ttm21f_vdwC=-2.147141323D+3,ttm21f_vdwD=1.d13, ttm21f_vdwE=13.2d0
double precision, parameter :: ttm21f_aDD=0.3d0, ttm21f_aCCaCD=0.2d0
double precision, parameter :: ttm21f_polfacO=0.837d0, ttm21f_polfacH = 0.496d0, ttm21f_polfacM=0.837d0
double precision, parameter :: ttm21f_polarO=0.837d0, ttm21f_polarH = 0.496d0, ttm21f_polarM=0.837d0
double precision, parameter :: ttm21f_gammaM=0.426706882d0
double precision, parameter :: ttm21f_dms_param1 = 0.d0, ttm21f_dms_param2 = 0.d0,  &
                               ttm21f_dms_param3 = 0.d0;

!----------------------------------------------------------------------------!
! Parameters of the TTM3F potential                                          !
!----------------------------------------------------------------------------!
double precision, parameter :: ttm3f_vdwA=0.d0,ttm3f_vdwB=0.d0,   &           
               ttm3f_vdwC=-0.72298855D+03,ttm3f_vdwD=0.10211829D+06, ttm3f_vdwE=0.37170376D+01 !(checked as Ang units) 
double precision, parameter :: ttm3f_aDD=0.175d0, ttm3f_aCCaCD=0.175d0
double precision, parameter :: ttm3f_polfacO=00.837d0, ttm3f_polfacH = 0.496d0, ttm3f_polfacM=0.837d0
double precision, parameter :: ttm3f_polarO=0.d0, ttm3f_polarH = 0.d0, ttm3f_polarM=1.444d0
double precision, parameter :: ttm3f_gammaM=0.46d0
double precision, parameter :: ttm3f_dms_param1 = 0.5d0, ttm3f_dms_param2 = 0.9578d0,  &
                               ttm3f_dms_param3 = 0.012d0;

!----------------------------------------------------------------------------!
! Parameters of the SPCfw potential                                          !
!----------------------------------------------------------------------------!
double precision, parameter ::    qspcfw_vdwA=629326.57249877361701443513d0, qspcfw_vdwB=0.d0, &
                       qspcfw_vdwC=-625.50206521141523463400d0, qspcfw_vdwD=0.d0, qspcfw_vdwE=0.d0, &
                        qspcfw_qO = -0.84d0, qspcfw_qH = 0.42d0, &
                        qspcfw_roh0=1.012d0,  qspcfw_th0=1.95476876223364912614d0, &  ! =112.0 deg
                        qspcfw_harmkb = 529.581d0, qspcfw_harmka = 37.95d0

!----------------------------------------------------------------------------!
! Parameters of the SPCf potential                                           !
!----------------------------------------------------------------------------!
double precision, parameter ::    spcf_vdwA=628213.72817108343399678768d0, spcf_vdwB=0.d0, &
                       spcf_vdwC=-624.97858519194071637260d0, spcf_vdwD=0.d0, spcf_vdwE=0.d0, &
                        spcf_qO = -0.82d0, spcf_qH = 0.41d0, &
                        spcf_roh0=1.d0,  spcf_th0=1.91061193215819258784d0, &  ! =112.0 deg (?)
                         spcf_Ap=671.288262912d0,  spcf_Bp=164.376d0, spcf_Cp=-211.536d0, spcf_Dp=111.744d0
!                        spcf_harmkb = 529.581d0, spcf_harmka = 37.95d0

double precision :: vdwA, vdwB, vdwC, vdwD, vdwE, aDD, accaCD, polarO, polarH, polarM, gammaM
double precision :: harmka, harmkb, roh0, th0, qO, qH

!----------------------------------------------------------------------------!
! Auxiliary variables and arrays                                             !
!----------------------------------------------------------------------------!
z

double precision, dimension(:,:), allocatable :: R, dR, dipt, dip, olddip, Efq, Efd
double precision, dimension(:), allocatable :: charge, phi
double precision, dimension(:,:,:,:), allocatable :: grdq
double precision, dimension(4) :: polfac
double precision, dimension(3,3) :: vir_rec
integer :: fO, lO, fH, lH, fM, lM, fO3, lO3, fH3, lH3, fM3, lM3, fd, ld, fd3, ld3
integer :: fdI, ldI
integer :: NatomsM
double precision :: Rc_sq, p4V, aewald
!double precision :: Uvdw_lrc0, Uvdw_lrc, Umon_tmp, Uelec_tmp, Uvdw_tmp, Uind_tmp
double precision :: Umon_tmp, Uelec_tmp, Uvdw_tmp, Uind_tmp, shiftA, shiftB
double precision :: shiftAc, shiftBc 

integer, dimension(:,:), allocatable :: k_hkl
double precision, dimension(:,:), allocatable :: coskx,cosky,coskz, sinkx,sinky,sinkz, sum_rec
double precision, dimension(:), allocatable   :: coskxky,coskr,   sinkxky,sinkr
double precision, dimension(:,:), allocatable :: kvec
double precision, dimension(:), allocatable :: rQq, iQq, rQm, iQm
double precision, dimension(:), allocatable :: Ak, dKk, dKkX, tmp_arr, dKkXY
integer :: nkvec, hmax, kmax, lmax, ewald_nk, ewald_ksqmax

double precision, dimension(:,:), allocatable :: fr, theta1, theta2, theta3,  dtheta1,dtheta2,dtheta3
integer, dimension(:,:), allocatable :: ifr
double precision, dimension(:,:), allocatable :: d2theta1,d2theta2,d2theta3
double precision, dimension(:), allocatable :: work
integer :: my_ndim3S, my_ndim3E, blk_ndim3 
integer :: my_ndim2S, my_ndim2E
logical :: polar_model


logical :: save_dd

end module pot_mod
