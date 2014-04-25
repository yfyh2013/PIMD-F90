module system_mod
   logical :: guess_initdip
   logical :: print_dipiters

   integer :: pot_model !  =2 for ttm21f and =3 for ttm3f

   double precision, dimension(:,:,:,:), allocatable :: previous_dips

   double precision, dimension(3) :: box, boxi
   double precision :: Umon, Uvdw, Uelec, Uind, Uvdw_lrc0, Uvdw_lrc

   double precision :: Rc, rc1, Rc2, volume, volume_init, eps_ewald
   double precision :: polar_sor, polar_eps
   integer :: polar_maxiter !changed to integer by D. Elton
   integer :: Natoms, Nwaters 
   integer, dimension(:,:), allocatable :: neigh_list
   integer :: predict_step
  double precision :: const_ts1, const_ts2, const_ts3

  double precision, dimension(:,:,:), allocatable :: tx_dip !moved here by Dan Elton

end module system_mod
