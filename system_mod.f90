module system_mod
   !storage space for global variables

   logical :: guess_initdip
   logical :: print_dipiters

   integer :: pot_model !  =2 for ttm21f and =3 for ttm3f

   double precision, dimension(3)   :: box, boxi
   double precision, dimension(3,3) :: siesta_box

   double precision :: Umon, Uvdw, Uelec, Uind, Uvdw_lrc0, Uvdw_lrc

   double precision :: Rc, rc1, Rc2, volume, volume_init, eps_ewald
   double precision :: polar_sor, polar_eps
   integer :: polar_maxiter !changed to integer by D. Elton
   integer :: Natoms, Nwaters, Nbeads
   integer, dimension(:,:), allocatable :: neigh_list
   integer :: predict_step
   double precision :: const_ts1, const_ts2, const_ts3

   logical :: CONTRACTION, MONOMERPIMD
   character(len=200) :: sys_label 

end module system_mod
