module m_data
  use m_cell
  implicit none
  integer :: testcase = 11
  !--------------------------------------------------------------------------
  ! Physical variables
  double precision :: elasticmodulus ! G>0 in viscoelasticity
  double precision :: oneoverell ! 1/ell>0 in viscoelasticity
  double precision :: oneoverlambda ! 1/lambda>0 in viscoelasticity
  double precision :: gavrilyuk ! 1 or 0: uses Gavrilyuk entrsophy variables or not
  double precision :: theta ! slope (gravity inclination in degrees: -90<theta<90)
  double precision :: Cfriction_u ! Navier friction coefficient
  double precision :: Cfriction_Phi ! microenstrophy friction coefficient in Gavrilyuk model
  !--------------------------------------------------------------------------
  integer :: I0 ! time-splitting
  integer :: I1 ! three-point numerical flux (order 1)
  integer :: I2 ! flux limiter (order 2 reconstruction) ; not used if 0
  double precision :: Tmax     ! final time
  double precision :: dtmin    ! minimal time-step tolerated
  integer :: Ntmax             ! maximal number of time iterations tolerated
  double precision, dimension(:,:), allocatable :: thist ! history of times and time steps
  integer :: Ntmax_clock       ! maximal number of times postprocessed: storage
  double precision :: dt_clock ! period of output/postprocessing dt<=dt_clock
  double precision :: t_clock  ! time of the next postprocessing
  integer :: Nx                ! number of spatial DOF (cells) and iterator
  double precision :: dx       ! cell volume (spatial step size)
  type (t_cell), dimension(:), allocatable :: cell
  double precision :: CFL      ! condition number
  double precision :: alphaspeed
end module m_data