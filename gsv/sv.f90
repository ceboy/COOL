program sv
  use m_zeromachine
  use m_physics
  use m_data
  implicit none
  !--------------------------------------------------------------------------
  integer :: Neq ! number of equations to be solved
  double precision :: dt    ! current time step
  integer :: nt = 0            ! time iteration number
  double precision :: t = 0.   ! current time
  integer :: ix                ! spatial DOF (cells) iterator
  type (t_cell), dimension(:), allocatable :: interf ! buffer cell
  double precision :: t_neighbour ! maximal time allowed by a CFL=1
  character (len=20) :: mystring ! buffer string
  logical :: mylogical
  integer :: myinteger, stencil
  !--------------------------------------------------------------------------
  interface
    subroutine riemann(t_neighbour,interf)
      use m_cell
      type (t_cell), dimension(:) :: interf 
      double precision :: t_neighbour
    end subroutine riemann
    subroutine printout(stencil,cell)
      use m_cell
      integer :: stencil
      type (t_cell), dimension(:) :: cell 
    end subroutine printout
    subroutine writeout(cell)
      use m_cell
      type (t_cell), dimension(:) :: cell 
    end subroutine writeout
  end interface
  ! Initialization ----------------------------------------------------------------
  myzeromachine = 1.
  do while(1./=1.+myzeromachine)
    myzeromachine = myzeromachine/2.
  end do
  print '("Zero machine = ",e12.6)', myzeromachine
  print *, 'Initialization'
  call svini ! <<<<<<<<<<<<<<<<<<<<<<<< INITIALIZE THE (GLOBAL) VARIABLES OF m_data
  print *, 'The step size in space is ', dx
  print *, 'Post-processing every ', dt_clock
  ! Parameterization ---------------------------------------------------------------
  open(unit=0,file='dt.res',form='formatted',status='new')
  open(unit=1,file='h.res',form='formatted',status='new')
  open(unit=2,file='Q.res',form='formatted',status='new')
  open(unit=3,file='u.res',form='formatted',status='new')
  open(unit=4,file='P.res',form='formatted',status='new')
  open(unit=5,file='phi.res',form='formatted',status='new')
  open(unit=11,file='sxx.res',form='formatted',status='new')
  open(unit=12,file='szz.res',form='formatted',status='new')
  ! >>> ABOVE: to be improved (writeout dependency on unit numbers)
  ! Post-processing and boundary conditions applied to copy interf -------------------
  allocate( interf(Nx+2) )
  interf = cell
  interf(1) = interf(2) ! no-flux
  interf(Nx+2) = interf(Nx+1) ! no-flux
  stencil = 5 ! (Nx+2)/2 ! number of cells shown on left and right boundaries <= (Nx+2)/2
  call printout(stencil,interf)
  print '("* Saving data at ",f16.8," >= ",f16.8)', t, t_clock
  call writeout(interf)
  t_clock = dt_clock ! next post-processing
  ! Computations -----------------------------------------------------------------
  print *, 'Entering time loop'
  time_loop : do while((t<Tmax).AND.(nt<Ntmax))
    nt = nt+1 
    ! Flux obtained from homogeneous (approximate) Riemann problems --------------
! celltrunc(interf) ! to handle vacuum in Riemann: set to zero non-sensible values ?
    call riemann(t_neighbour,interf) !! Suliciu : choice ?
! interf = fluxtrunc(interf) ! to handle vacuum in fluxes: set to zero non-sensible values ?
    ! >>>> ABOVE: one should better create a new type "flux" with same attributes as cell !
    ! Time-step ------------------------------------------------------------------
    dt = t_neighbour*CFL
    if (t_neighbour<dtmin) then
      dt = dtmin
      print *, 'dt too small'
      exit
    end if
    mylogical = .FALSE.
    if ((t+dt>Tmax).OR.(t+dt>t_clock)) then
      mylogical = .TRUE.
      dt = min(Tmax-t,t_clock-t) !dt = minval( (\ Tmax-t, dt_clock \) ) 
    end if
    print '("Time step at iteration ",i10," = ",e10.3," time = ",f10.2)',nt,dt,t
    thist(1,nt) = t ; thist(2,nt) = t_neighbour*CFL ; thist(3,nt) = dt
    write( 0, '(3e15.6)') thist(1:3,nt)
    t = t+dt
    ! First-order time-splitting ------------------------------------------------
    cell(2:Nx+1)%depth = cell(2:Nx+1)%depth + dt*interf(2:Nx+1)%depth
    if(minval(cell(2:Nx+1)%depth)<0.) then
      print *, 'Negative depth!'; exit
    end if
    cell(2:Nx+1)%discharge = cell(2:Nx+1)%discharge + dt*interf(2:Nx+1)%discharge
    cell(2:Nx+1)%htracer = cell(2:Nx+1)%htracer + dt*interf(2:Nx+1)%htracer
    cell(2:Nx+1)%hsigmaxx = cell(2:Nx+1)%hsigmaxx + dt*interf(2:Nx+1)%hsigmaxx
    cell(2:Nx+1)%hsigmazz = cell(2:Nx+1)%hsigmazz + dt*interf(2:Nx+1)%hsigmazz
    !where( (cell(2:Nx+1)%depth>0.) ) 
    where( (cell(2:Nx+1)%depth>myzeromachine*maxval(cell(2:Nx+1)%depth)) )
      cell(2:Nx+1)%velocity = cell(2:Nx+1)%discharge/cell(2:Nx+1)%depth
      cell(2:Nx+1)%tracer = cell(2:Nx+1)%htracer/cell(2:Nx+1)%depth
      cell(2:Nx+1)%sigmaxx = cell(2:Nx+1)%hsigmaxx/cell(2:Nx+1)%depth
      cell(2:Nx+1)%sigmazz = cell(2:Nx+1)%hsigmazz/cell(2:Nx+1)%depth
      cell(2:Nx+1)%pressure = g*cell(2:Nx+1)%depth**2/2 &
        + elasticmodulus*(cell(2:Nx+1)%sigmazz-cell(2:Nx+1)%sigmaxx)
    elsewhere
      cell(2:Nx+1)%velocity = 0.
      cell(2:Nx+1)%tracer = 0.
      cell(2:Nx+1)%sigmaxx = 0.
      cell(2:Nx+1)%sigmazz = 0.
      cell(2:Nx+1)%pressure = 0.
    end where
    ! Post-processing and boundary conditions applied to copy interf --------------
    interf = cell 
    interf(1) = interf(2) ! no-flux
    interf(Nx+2) = interf(Nx+1) ! no-flux
    if(mylogical) then
      call printout(stencil,interf)
      print '("* Saving data at ",f16.8," >= ",f16.8)', t, t_clock
      call writeout(interf)
      t_clock = t_clock + dt_clock
    end if
  end do time_loop 
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  function celltrunc(cell)
    use m_cell
    use m_zeromachine
    implicit none
    type (t_cell), dimension(:) :: cell
    type (t_cell), dimension(:), allocatable :: celltrunc
    double precision :: ref
    allocate( celltrunc(size(cell)) )
    celltrunc = cell
    where( celltrunc%depth<maxval(celltrunc%depth)*myzeromachine )
      celltrunc%depth = 0.
      celltrunc%discharge = 0.
      celltrunc%htracer = 0.
      celltrunc%hsigmaxx = 0.
      celltrunc%hsigmazz = 0.
      celltrunc%velocity = 0.
      celltrunc%tracer = 0.
      celltrunc%sigmaxx = 0.
      celltrunc%sigmazz = 0.
      celltrunc%pressure = 0.
    end where
  endfunction celltrunc
  function fluxtrunc(interf)
    use m_cell
    use m_zeromachine
    implicit none
    type (t_cell), dimension(:) :: interf
    type (t_cell), dimension(:), allocatable :: fluxtrunc
    double precision :: ref
    allocate( fluxtrunc(size(interf)) )
    fluxtrunc = interf
    ref = max(maxval(interf%depth),-minval(interf%depth))
    where( (interf%depth<=myzeromachine*ref).AND.(interf%depth>0.) ) fluxtrunc%depth=0.
    where( (interf%depth>=-myzeromachine*ref).AND.(interf%depth<0.) ) fluxtrunc%depth=0.
    ref = max(maxval(interf%discharge),-minval(interf%discharge))
    where( (interf%discharge<=myzeromachine*ref).AND.(interf%discharge>0.) ) fluxtrunc%discharge=0.
    where( (interf%discharge>=-myzeromachine*ref).AND.(interf%discharge<0.) ) fluxtrunc%discharge=0.
    ref = max(maxval(interf%htracer),-minval(interf%htracer))
    where( (interf%htracer<=myzeromachine*ref).AND.(interf%htracer>0.) ) fluxtrunc%htracer=0.
    where( (interf%htracer>=-myzeromachine*ref).AND.(interf%htracer<0.) ) fluxtrunc%htracer=0.
    ref = max(maxval(interf%hsigmaxx),-minval(interf%hsigmaxx))
    where( (interf%hsigmaxx<=myzeromachine*ref).AND.(interf%hsigmaxx>0.) ) fluxtrunc%hsigmaxx=0.
    where( (interf%hsigmaxx>=-myzeromachine*ref).AND.(interf%hsigmaxx<0.) ) fluxtrunc%hsigmaxx=0.
    ref = max(maxval(interf%hsigmazz),-minval(interf%hsigmazz))
    where( (interf%hsigmazz<=myzeromachine*ref).AND.(interf%hsigmazz>0.) ) fluxtrunc%hsigmazz=0.
    where( (interf%hsigmazz>=-myzeromachine*ref).AND.(interf%hsigmazz<0.) ) fluxtrunc%hsigmazz=0.
    return
  end function
end program sv

! ------------------------------------------------------------

subroutine riemann(t_neighbour,cell)
  use m_physics
  use m_cell
  implicit none
  type (t_cell), dimension(:), intent(inout) :: cell 
  double precision, intent(out) :: t_neighbour
  integer :: N0, Nface
  double precision, dimension(:,:), allocatable :: fluxleft, fluxright
  !!! * Suliciu relaxation solver: diagonal system with additional variables
  double precision, dimension(:), allocatable :: leftspeed, rightspeed  
  double precision, dimension(:), allocatable :: leftparameter, rightparameter  
  double precision, dimension(:), allocatable :: leftdepth, rightdepth
  double precision, dimension(:), allocatable :: suliciuvelocity, suliciupressure
  double precision, dimension(:), allocatable :: leftvelocity, rightvelocity
!   double precision, dimension(:), allocatable :: leftsigmaxx, rightsigmaxx
!   double precision, dimension(:), allocatable :: leftsigmazz, rightsigmazz
  double precision, dimension(:), allocatable :: lefthsigmaxx, righthsigmaxx
  double precision, dimension(:), allocatable :: lefthsigmazz, righthsigmazz
  logical, dimension(:), allocatable :: mylogical
  integer :: ii = 0
  !--------------------------------------------------------------------------
  ! interface
  !   subroutine suliciu_initialization( leftparameter, rightparameter, cell ) 
  !     use m_cell
  !     double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
  !     type(t_cell), dimension(:), intent(in) :: cell
  !   end subroutine suliciu_initialization
  ! end interface
  !--------------------------------------------------------------------------
  N0 = size(cell)
  Nface = N0-1
  allocate( fluxleft(5,Nface), fluxright(5,Nface), &
    leftspeed(Nface), rightspeed(Nface), &
    leftparameter(Nface), rightparameter(Nface), &
    leftdepth(Nface), rightdepth(Nface), &
    suliciuvelocity(Nface), suliciupressure(Nface), &
    leftvelocity(Nface), rightvelocity(Nface), &
    mylogical(Nface) &
  )
  !!! * explicit initialization of the parameter
  leftspeed = sqrt(g*cell(1:N0-1)%depth+elasticmodulus*(3*cell(1:N0-1)%sigmazz+cell(1:N0-1)%sigmaxx))
  rightspeed = sqrt(g*cell(2:N0)%depth+elasticmodulus*(3*cell(2:N0)%sigmazz+cell(2:N0)%sigmaxx))
  call suliciu_initialization( leftspeed, rightspeed, cell )
  !!! * extremal wave speeds of Suliciu system
  leftvelocity = cell(1:N0-1)%velocity - leftspeed ! ==suliciuvelocity-leftspeed*(cell(1:N0-1)%depth/leftdepth)
  rightvelocity = cell(2:N0)%velocity + rightspeed ! ==suliciuvelocity+rightspeed*(cell(2:N0)%depth/rightdepth)
  !!! * intermediate variables
  leftparameter = leftspeed*cell(1:N0-1)%depth ; rightparameter = rightspeed*cell(2:N0)%depth
  !!! * Suliciu explicit solution
  suliciuvelocity = 0.
  suliciupressure = 0.
  leftdepth = 0.
  rightdepth = 0.
!   allocate( leftsigmaxx(Nface), rightsigmaxx(Nface), leftsigmazz(Nface), rightsigmazz(Nface) )
!   leftsigmaxx = 0.
!   leftsigmazz = 0.
!   rightsigmaxx = 0.
!   rightsigmazz = 0.
!! ABOVE: should we use intensive variables (e.g. to take into account some truncation ?) ?????
  allocate( lefthsigmaxx(Nface), righthsigmaxx(Nface), lefthsigmazz(Nface), righthsigmazz(Nface) )
  lefthsigmaxx = 0.
  lefthsigmazz = 0.
  righthsigmaxx = 0.
  righthsigmazz = 0.
  where( (cell(1:N0-1)%depth/=0.).AND.(cell(2:N0)%depth==0.) ) 
    suliciuvelocity = cell(1:N0-1)%velocity + cell(1:N0-1)%pressure/leftparameter
    ! suliciupressure = 0.
    leftdepth = cell(1:N0-1)%depth/( 1. + cell(1:N0-1)%pressure/(leftspeed*leftparameter) )
    ! rightdepth = 0.
    lefthsigmaxx = cell(1:N0-1)%hsigmaxx*(1. + cell(1:N0-1)%pressure/(leftspeed*leftparameter))
    lefthsigmazz = cell(1:N0-1)%hsigmazz/(1. + cell(1:N0-1)%pressure/(leftspeed*leftparameter))**3
  end where
  where( (cell(1:N0-1)%depth==0.).AND.(cell(2:N0)%depth/=0.) )
    suliciuvelocity = cell(2:N0)%velocity - cell(2:N0)%pressure/rightparameter
    ! suliciupressure = 0.
    ! leftdepth = 0.
    rightdepth = cell(2:N0)%depth/( 1. + cell(2:N0)%pressure/(rightspeed*rightparameter) )
    righthsigmaxx = cell(2:N0)%hsigmaxx*(1. + cell(2:N0)%pressure/(rightspeed*rightparameter))
    righthsigmazz = cell(2:N0)%hsigmazz/(1. + cell(2:N0)%pressure/(rightspeed*rightparameter))**3
  end where
  where( (cell(1:N0-1)%depth/=0.).AND.(cell(2:N0)%depth/=0.) )
!     suliciuvelocity = ( leftspeed*cell(1:N0-1)%discharge + rightspeed*cell(2:N0)%discharge &
!                         + cell(1:N0-1)%pressure - cell(2:N0)%pressure ) &
!                      /( leftspeed*cell(1:N0-1)%depth + rightspeed*cell(2:N0)%depth )
    suliciuvelocity = ( leftspeed*cell(1:N0-1)%discharge + rightspeed*cell(2:N0)%discharge &
                        + cell(1:N0-1)%pressure - cell(2:N0)%pressure ) &
                     /( rightparameter+leftparameter )
!     suliciuvelocity = ( leftparameter*cell(1:N0-1)%velocity + rightparameter*cell(2:N0)%velocity &
!       + cell(1:N0-1)%pressure - cell(2:N0)%pressure )/(rightparameter+leftparameter)
!! ABOVE: should we use velocity (e.g. to take into account some truncation ?) ?????
!     suliciupressure = ( cell(2:N0)%pressure*leftspeed*cell(1:N0-1)%depth &
!                         + cell(1:N0-1)%pressure*rightspeed*cell(2:N0)%depth &
!     + leftspeed*rightspeed*(cell(1:N0-1)%discharge*cell(2:N0)%depth-cell(2:N0)%discharge*cell(1:N0-1)%depth) ) &
!                      /( leftspeed*cell(1:N0-1)%depth + rightspeed*cell(2:N0)%depth )
    suliciupressure = ( cell(2:N0)%pressure*leftparameter + cell(1:N0-1)%pressure*rightparameter &
     + leftspeed*rightspeed*(cell(1:N0-1)%discharge*cell(2:N0)%depth-cell(2:N0)%discharge*cell(1:N0-1)%depth) ) &
                     /( rightparameter+leftparameter )
!     suliciupressure = ( cell(2:N0)%pressure*leftparameter + cell(1:N0-1)%pressure*rightparameter &
!       + leftparameter*rightparameter*(cell(1:N0-1)%velocity-cell(2:N0)%velocity) )/(rightparameter+leftparameter)
!! ABOVE: should we use velocity (e.g. to take into account some truncation ?) ?????
    leftdepth = cell(1:N0-1)%depth/(1.-((suliciupressure-cell(1:N0-1)%pressure)/(leftspeed*leftparameter)))
    rightdepth = cell(2:N0)%depth/(1.-((suliciupressure-cell(2:N0)%pressure)/(rightspeed*rightparameter)))
!     leftsigmaxx = cell(1:N0-1)%sigmaxx*(cell(1:N0-1)%depth/leftdepth)**2
!     leftsigmazz = cell(1:N0-1)%sigmazz*(leftdepth/cell(1:N0-1)%depth)**2
!     lefthsigmaxx = cell(1:N0-1)%hsigmaxx*(cell(1:N0-1)%depth/leftdepth)
!     lefthsigmazz = cell(1:N0-1)%hsigmazz*(leftdepth/cell(1:N0-1)%depth)**3
    lefthsigmaxx = cell(1:N0-1)%hsigmaxx*(1.-((suliciupressure-cell(1:N0-1)%pressure)/(leftspeed*leftparameter)))
    lefthsigmazz = cell(1:N0-1)%hsigmazz/(1.-((suliciupressure-cell(1:N0-1)%pressure)/(leftspeed*leftparameter)))**3
!     rightsigmaxx = cell(2:N0)%sigmaxx*(cell(2:N0)%depth/rightdepth)**2
!     rightsigmazz = cell(2:N0)%sigmazz*(rightdepth/cell(2:N0)%depth)**2
!     righthsigmaxx = cell(2:N0)%hsigmaxx*(cell(2:N0)%depth/rightdepth)
!     righthsigmazz = cell(2:N0)%hsigmazz*(rightdepth/cell(2:N0)%depth)**3
    righthsigmaxx = cell(2:N0)%hsigmaxx*(1.-((suliciupressure-cell(2:N0)%pressure)/(rightspeed*rightparameter)))
    righthsigmazz = cell(2:N0)%hsigmazz/(1.-((suliciupressure-cell(2:N0)%pressure)/(rightspeed*rightparameter)))**3
  end where
!   print *, '---------DEBUG-------------'
!   print *, N0
!   mylogical = .FALSE.
!   where( suliciuvelocity>rightvelocity )
!     mylogical = .TRUE.
!   end where
!   where( suliciuvelocity<leftvelocity )
!     mylogical = .TRUE.
!   end where
!   mylogical = (cell(1:N0-1)%depth==0.).AND.(cell(2:N0)%depth==0.)
!   ii = 1
!   do while( mylogical(ii).eqv..FALSE. )
!     ii = ii+1
!   end do
!   print *, ii
!   print *, cell(ii-1)%depth
!   print *, cell(ii)%depth
!   print *, leftspeed
!   print *, rightspeed
!   print *, leftparameter
!   print *, rightparameter
!   print '(21f10.3)', leftvelocity
!   print '(20f10.3)', cell%velocity
!   print '(21f10.3)', suliciuvelocity
!   print '(21f10.3)', rightvelocity
!   print '(20f10.3)', cell%pressure
!   print '(21f10.3)', suliciupressure
!   print *, cell%depth
!   print *, leftdepth
!   print *, rightdepth
  print *, cell%hsigmaxx*cell%depth
  print *, lefthsigmaxx*leftdepth
  print *, righthsigmaxx*rightdepth
  print *, cell%hsigmazz/(1.e-16+cell%depth)**3
  print *, lefthsigmazz/(1.e-16+leftdepth)**3
  print *, righthsigmazz/(1.e-16+rightdepth)**3
!   print *, '---------DEBUG-------------'
  !!! * Suliciu flux computation
  where( (suliciuvelocity>=0.).AND.(leftvelocity>=0.) )
    fluxleft(1,:) = cell(1:N0-1)%depth*cell(1:N0-1)%velocity
    fluxleft(2,:) = cell(1:N0-1)%depth*cell(1:N0-1)%velocity**2 + cell(1:N0-1)%pressure
!     fluxleft(3,:) = cell(1:N0-1)%tracer*cell(1:N0-1)%discharge
!     fluxleft(4,:) = cell(1:N0-1)%sigmaxx*cell(1:N0-1)%discharge
!     fluxleft(5,:) = cell(1:N0-1)%sigmazz*cell(1:N0-1)%discharge
    fluxleft(3,:) = cell(1:N0-1)%htracer*cell(1:N0-1)%velocity
    fluxleft(4,:) = cell(1:N0-1)%hsigmaxx*cell(1:N0-1)%velocity
    fluxleft(5,:) = cell(1:N0-1)%hsigmazz*cell(1:N0-1)%velocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = cell(2:N0)%sigmaxx*cell(2:N0)%discharge & 
!       - leftvelocity*(leftsigmaxx*leftdepth-cell(1:N0-1)%sigmaxx*cell(1:N0-1)%depth) &
!       - suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx*leftdepth) &
!       - rightvelocity*(cell(2:N0)%sigmaxx*cell(2:N0)%depth-rightsigmaxx*rightdepth)
!     fluxright(5,:) = cell(2:N0)%sigmazz*cell(2:N0)%discharge & 
!       - leftvelocity*(leftsigmazz*leftdepth-cell(1:N0-1)%sigmazz*cell(1:N0-1)%depth) &
!       - suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth) &
!       - rightvelocity*(cell(2:N0)%sigmazz*cell(2:N0)%depth-rightsigmazz*rightdepth)
    fluxright(4,:) = cell(2:N0)%hsigmaxx*cell(2:N0)%velocity & 
      - rightvelocity*(cell(2:N0)%hsigmaxx-righthsigmaxx) &
      - suliciuvelocity*(righthsigmaxx-lefthsigmaxx) &
      - leftvelocity*(lefthsigmaxx-cell(1:N0-1)%hsigmaxx)
    fluxright(5,:) = cell(2:N0)%hsigmazz*cell(2:N0)%velocity & 
      - rightvelocity*(cell(2:N0)%hsigmazz-righthsigmazz) &
      - suliciuvelocity*(righthsigmazz-lefthsigmazz) &
      - leftvelocity*(lefthsigmazz-cell(1:N0-1)%hsigmazz)
  end where
  where( (suliciuvelocity>=0.).AND.(leftvelocity<0.) )
    fluxleft(1,:) = leftdepth*suliciuvelocity
    fluxleft(2,:) = leftdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = cell(1:N0-1)%tracer*leftdepth*suliciuvelocity
!     fluxleft(4,:) = cell(1:N0-1)%sigmaxx*cell(1:N0-1)%discharge &
!       + leftvelocity*(leftsigmaxx*leftdepth-cell(1:N0-1)%sigmaxx*cell(1:N0-1)%depth)
!     fluxleft(5,:) = cell(1:N0-1)%sigmazz*cell(1:N0-1)%discharge &
!       + leftvelocity*(leftsigmazz*leftdepth-cell(1:N0-1)%sigmazz*cell(1:N0-1)%depth)
    fluxleft(4,:) = cell(1:N0-1)%hsigmaxx*cell(1:N0-1)%velocity &
      + leftvelocity*(lefthsigmaxx-cell(1:N0-1)%hsigmaxx)
    fluxleft(5,:) = cell(1:N0-1)%hsigmazz*cell(1:N0-1)%velocity &
      + leftvelocity*(lefthsigmazz-cell(1:N0-1)%hsigmazz)
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = cell(2:N0)%sigmaxx*cell(2:N0)%discharge & 
!       - suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx*leftdepth) &
!       - rightvelocity*(cell(2:N0)%sigmaxx*cell(2:N0)%depth-rightsigmaxx*rightdepth)
!     fluxright(5,:) = cell(2:N0)%sigmazz*cell(2:N0)%discharge & 
!       - suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth) &
!       - rightvelocity*(cell(2:N0)%sigmazz*cell(2:N0)%depth-rightsigmazz*rightdepth)
    fluxright(4,:) = cell(2:N0)%hsigmaxx*cell(2:N0)%velocity & 
      - rightvelocity*(cell(2:N0)%hsigmaxx-righthsigmaxx) &
      - suliciuvelocity*(righthsigmaxx-lefthsigmaxx)
    fluxright(5,:) = cell(2:N0)%hsigmazz*cell(2:N0)%velocity & 
      - rightvelocity*(cell(2:N0)%hsigmazz-righthsigmazz) &
      - suliciuvelocity*(righthsigmazz-lefthsigmazz)
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity>=0.) )
    fluxleft(1,:) = rightdepth*suliciuvelocity
    fluxleft(2,:) = rightdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = cell(2:N0)%tracer*rightdepth*suliciuvelocity
!     fluxleft(4,:) = cell(1:N0-1)%sigmaxx*cell(1:N0-1)%discharge &
!       + leftvelocity*(leftsigmaxx*leftdepth-cell(1:N0-1)%sigmaxx*cell(1:N0-1)%depth) &
!       + suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx)
!     fluxleft(5,:) = cell(1:N0-1)%sigmazz*cell(1:N0-1)%discharge &
!       + leftvelocity*(leftsigmazz*leftdepth-cell(1:N0-1)%sigmazz*cell(1:N0-1)%depth) &
!       + suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth)
    fluxleft(4,:) = cell(1:N0-1)%hsigmaxx*cell(1:N0-1)%velocity &
      + leftvelocity*(lefthsigmaxx-cell(1:N0-1)%hsigmaxx) &
      + suliciuvelocity*(righthsigmaxx-lefthsigmaxx)
    fluxleft(5,:) = cell(1:N0-1)%hsigmazz*cell(1:N0-1)%velocity &
      + leftvelocity*(lefthsigmazz-cell(1:N0-1)%hsigmazz) &
      + suliciuvelocity*(righthsigmazz-lefthsigmazz)
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = cell(2:N0)%sigmaxx*cell(2:N0)%discharge & 
!       - rightvelocity*(cell(2:N0)%sigmaxx*cell(2:N0)%depth-rightsigmaxx*rightdepth)
!     fluxright(5,:) = cell(2:N0)%sigmazz*cell(2:N0)%discharge & 
!       - rightvelocity*(cell(2:N0)%sigmazz*cell(2:N0)%depth-rightsigmazz*rightdepth)
    fluxright(4,:) = cell(2:N0)%hsigmaxx*cell(2:N0)%velocity & 
      - rightvelocity*(cell(2:N0)%hsigmaxx-righthsigmaxx)
    fluxright(5,:) = cell(2:N0)%hsigmazz*cell(2:N0)%velocity & 
      - rightvelocity*(cell(2:N0)%hsigmazz-righthsigmazz)
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity<0.) )
    fluxleft(1,:) = cell(2:N0)%depth*cell(2:N0)%velocity
    fluxleft(2,:) = cell(2:N0)%depth*cell(2:N0)%velocity**2 + cell(2:N0)%pressure 
!     fluxleft(3,:) = cell(2:N0)%tracer*cell(2:N0)%discharge
!     fluxleft(4,:) = cell(1:N0-1)%sigmaxx*cell(1:N0-1)%discharge &
!       + leftvelocity*(leftsigmaxx*leftdepth-cell(1:N0-1)%sigmaxx*cell(1:N0-1)%depth) &
!       + suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx*leftdepth) &
!       + rightvelocity*(cell(2:N0)%sigmaxx*cell(2:N0)%depth-rightsigmaxx*rightdepth)
!     fluxleft(5,:) = cell(1:N0-1)%sigmazz*cell(1:N0-1)%discharge &
!       + leftvelocity*(leftsigmazz*leftdepth-cell(1:N0-1)%sigmazz*cell(1:N0-1)%depth) &
!       + suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth) &
!       + rightvelocity*(cell(2:N0)%sigmazz*cell(2:N0)%depth-rightsigmazz*rightdepth)
    fluxleft(3,:) = cell(2:N0)%htracer*cell(2:N0)%velocity
    fluxleft(4,:) = cell(1:N0-1)%hsigmaxx*cell(1:N0-1)%velocity &
      + leftvelocity*(lefthsigmaxx-cell(1:N0-1)%hsigmaxx) &
      + suliciuvelocity*(righthsigmaxx-lefthsigmaxx) &
      + rightvelocity*(cell(2:N0)%hsigmaxx-righthsigmaxx)
    fluxleft(5,:) = cell(1:N0-1)%hsigmazz*cell(1:N0-1)%velocity &
      + leftvelocity*(lefthsigmazz-cell(1:N0-1)%hsigmazz) &
      + suliciuvelocity*(righthsigmazz-lefthsigmazz) &
      + rightvelocity*(cell(2:N0)%hsigmazz-righthsigmazz)
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = cell(2:N0)%sigmaxx*cell(2:N0)%discharge
!     fluxright(5,:) = cell(2:N0)%sigmazz*cell(2:N0)%discharge
    fluxright(4,:) = cell(2:N0)%hsigmaxx*cell(2:N0)%velocity
    fluxright(5,:) = cell(2:N0)%hsigmazz*cell(2:N0)%velocity
  end where
!   print *, '---------DEBUG-------------'
!   print '(21f10.3)', fluxleft(1,:)
!   print '(21f10.3)', fluxleft(2,:)
!   print '(21f10.3)', fluxleft(3,:)
!   print '(21f10.3)', fluxleft(4,:)
!   print '(21f10.3)', fluxleft(5,:)
!   print *, '---------DEBUG-------------'
  cell(2:N0-1)%depth = (fluxright(1,1:Nface-1)-fluxleft(1,2:Nface))/cell(2:N0-1)%volume
  cell(2:N0-1)%discharge = (fluxright(2,1:Nface-1)-fluxleft(2,2:Nface))/cell(2:N0-1)%volume
  cell(2:N0-1)%htracer = (fluxright(3,1:Nface-1)-fluxleft(3,2:Nface))/cell(2:N0-1)%volume
  cell(2:N0-1)%hsigmaxx = (fluxright(4,1:Nface-1)-fluxleft(4,2:Nface))/cell(2:N0-1)%volume
  cell(2:N0-1)%hsigmazz = (fluxright(5,1:Nface-1)-fluxleft(5,2:Nface))/cell(2:N0-1)%volume
  where( rightvelocity<0 ) rightvelocity = -rightvelocity
  where( leftvelocity<0 ) leftvelocity = -leftvelocity
  t_neighbour = minval( (/cell(2:N0)%volume/rightvelocity,cell(1:N0-1)%volume/leftvelocity/) )
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  subroutine suliciu_initialization( leftparameter, rightparameter, cell ) 
    use m_cell
    implicit none
    double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
    type(t_cell), dimension(:), intent(in) :: cell
    integer :: N0, Nface
    double precision, dimension(:), allocatable :: temp0, temp1, temp2, temp3
    double precision :: alpha
    N0 = size(cell)
    Nface = N0-1
    allocate( temp0(Nface), temp1(Nface), temp2(Nface), temp3(Nface) )
    temp0 = rightparameter*cell(2:N0)%depth + leftparameter*cell(1:N0-1)%depth
    temp1 = cell(2:N0)%pressure-cell(1:N0-1)%pressure ! Pr-Pl
    where(temp0==0.)
      temp2 = 0.
      temp3 = 0.
    end where
    where((temp1>=0.).AND.(temp0/=0.))
      temp2 = temp1/temp0
      temp3 = 0.
    end where
    where((temp1<0.).AND.(temp0/=0.))
      temp2 = 0.
      temp3 = -temp1/temp0
    end where
    temp1 = cell(1:N0-1)%velocity-cell(2:N0)%velocity ! u_l-u_r
    where(temp1<0.) temp1=0.
    temp2 = ( temp2 + temp1 )
    temp3 = ( temp3 + temp1 )
    ! --
    alpha = 2. !1.5
    !leftparameter = leftparameter*(1. + alpha*temp2/leftparameter)
    !rightparameter = rightparameter*(1. + alpha*temp3/rightparameter)
    leftparameter = leftparameter + alpha*temp2
    rightparameter = rightparameter + alpha*temp3
  end subroutine suliciu_initialization
end subroutine riemann

! ------------------------------------------------------------

subroutine printout(stencil,cell)
  use m_cell
  implicit none
  integer :: stencil, ix, Nx
  type (t_cell), dimension(:) :: cell 
  character (len=20) :: mystring ! buffer string
  Nx = size(cell)-2
  print *, 'Cell values including boundary ghost cells'
  write(mystring,*) stencil
  print '('//mystring//'f16.8,"   ...",'//mystring//'f16.8)', &
    (/ (cell(ix)%depth, ix=1,stencil) , &
       (cell(ix)%depth, ix=Nx+2-stencil+1,Nx+2) /) 
!   print '('//mystring//'f16.8,"   ...",'//mystring//'f16.8)', &
!     (/ (cell(ix)%tracer, ix=1,stencil) , &
!         (cell(ix)%tracer, ix=Nx+2-stencil+1,Nx+2) /) 
!   print '('//mystring//'f16.8,"   ...",'//mystring//'f16.8)', &
!     (/ (cell(ix)%tracer/cell(ix)%depth, ix=1,stencil) , &
!        (cell(ix)%tracer/cell(ix)%depth, ix=Nx+2-stencil+1,Nx+2) /) 
  print '('//mystring//'f16.8,"   ...",'//mystring//'f16.8)', &
    (/ (cell(ix)%sigmazz, ix=1,stencil) , &
       (cell(ix)%sigmazz, ix=Nx+2-stencil+1,Nx+2) /) 
end subroutine printout

subroutine writeout(cell)
  use m_cell
  implicit none
  type (t_cell), dimension(:) :: cell 
  integer :: Nx
  character (len=20) :: mystring ! buffer string
  Nx = size(cell)-2
  write(mystring,*) Nx
  write( 1, '('//mystring//'f16.8)') cell(2:(Nx+1))%depth
  write( 2, '('//mystring//'f16.8)') cell(2:(Nx+1))%discharge
  write( 3, '('//mystring//'f16.8)') cell(2:(Nx+1))%velocity
  write( 5, '('//mystring//'f16.8)') cell(2:(Nx+1))%tracer
  write( 11, '('//mystring//'f16.8)') cell(2:(Nx+1))%sigmaxx
  write( 12, '('//mystring//'f16.8)') cell(2:(Nx+1))%sigmazz
  write( 4, '('//mystring//'f16.8)') cell(2:(Nx+1))%pressure
end subroutine writeout


