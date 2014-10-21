! module pour post processing ??
program sv
  use m_data
  implicit none
  !--------------------------------------------------------------------------
  double precision :: dt    ! current time step
  integer :: nt = 0            ! time iteration number
  double precision :: t = 0.   ! current time
  integer :: ix                ! spatial DOF (cells) iterator
  type (t_cell), dimension(:), allocatable :: cell1 ! buffer cell
  integer, dimension(:), allocatable :: indicesx
  character (len=20) :: mystring ! buffer string
  logical :: mylogical
  integer :: myinteger
  !--------------------------------------------------------------------------
  interface
    subroutine riemann(cell0,CFL,dtmin,dtmax,dt,g)
      use m_cell
      double precision :: CFL,dtmin,dtmax,g
      double precision :: dt
      type (t_cell), dimension(:) :: cell0 
    end subroutine
    subroutine printout(stencil)
      use m_data
      integer :: stencil
      character (len=20) :: mystring
    end subroutine printout
  end interface
  ! Initialization ----------------------------------------------------------------
  print *, 'Initialization'
  ! print *, 'Number of equidistant points...?'
  ! read *, Nx ! surcharge svini() ??
  call svini ! <<<<<<<<<<<<<<<<<<<<<<<< INITIALIZE THE (GLOBAL) VARIABLES OF m_data
  print *, 'The step size in space is ', dx
  print *, 'Post-processing every ', dt_clock
  allocate( indicesx(Nx+2) )
  ! Parameterization ---------------------------------------------------------------
  open(unit=1,file='h.res',form='formatted',status='new')
  open(unit=2,file='Q.res',form='formatted',status='new')
  open(unit=3,file='u.res',form='formatted',status='new')
  open(unit=4,file='P.res',form='formatted',status='new')
  ! Post-processing -------------------------------------------------------------
  call printout(4)
  print '("* Saving data at ",f6.2)', t
  write(mystring, *) Nx
  write( 1, '('//mystring//'f6.2)') cell(2:(Nx+1))%depth
  write( 2, '('//mystring//'f6.2)') cell(2:(Nx+1))%discharge
  write( 3, '('//mystring//'f6.2)') cell(2:(Nx+1))%velocity
  write( 4, '('//mystring//'f6.2)') cell(2:(Nx+1))%pressure
  t_clock = dt_clock
  ! Computations ----------------------------------------------------------------
  print *, 'Entering time loop'
  allocate( cell1(Nx+2) )
  time_loop : do while(t<Tmax)
    nt = nt+1 
    if(nt==Ntmax) exit
    ! Boundary Conditions applied to copy cell1 ---------------------------------
    cell1 = cell
    cell1(1) = cell1(2) ! no-flux
    cell1(Nx+2) = cell1(Nx+1) ! no-flux
    ! Homogeneous Riemann problem solved (approximately) ------------------------
    call riemann(cell1,CFL,dtmin,min(Tmax-t,dt_clock),dt,g) !! Suliciu relaxation : add choice ?
    mylogical = (dt==min(Tmax-t,dt_clock))
    ! First-order time-splitting ------------------------------------------------
    cell = cell1
    ! Post-processing -----------------------------------------------------------
    print '("Time step at iteration",i10," = ",e10.3," time = ",f10.2)',nt,dt,t
    thist(1,nt) = t
    thist(2,nt) = dt
    t = t+dt
    if((t.ge.t_clock).OR.mylogical) then
      call printout(4)
      print '("* Saving data at ",f6.2," >= ",f6.2)', t, t_clock
      write(mystring,*) Nx
      write( 1, '('//mystring//'f6.2)') cell(2:(Nx+1))%depth
      write( 2, '('//mystring//'f6.2)') cell(2:(Nx+1))%discharge
      write( 3, '('//mystring//'f6.2)') cell(2:(Nx+1))%velocity
      write( 4, '('//mystring//'f6.2)') cell(2:(Nx+1))%pressure
      t_clock = t_clock + dt_clock
    end if
  end do time_loop 
end program sv

! ------------------------------------------------------------

subroutine riemann(cell0,CFL,dtmin,dtmax,dt,g) !! Suliciu relaxation : add choice ?
  use m_cell
  implicit none
  double precision, intent(in) :: CFL,dtmin,dtmax,g
  double precision, intent(inout) :: dt
  type (t_cell), dimension(:), intent(inout) :: cell0 
  integer :: N0
  integer :: iface, Nface
  ! Suliciu relaxation solver: diagonal system with additional variables
  double precision, dimension(:), allocatable :: leftpressure, rightpressure, deltap
  double precision, dimension(:), allocatable :: leftparameter, rightparameter  
  double precision, dimension(:), allocatable :: leftdepth, rightdepth
  double precision, dimension(:), allocatable :: suliciuvelocity, suliciupressure
  double precision, dimension(:), allocatable :: leftvelocity, rightvelocity
  double precision, dimension(:,:), allocatable :: fluxleft, fluxright  
character*20 :: mystring
  N0 = size(cell0)
  Nface = N0-1
write(mystring, *) Nface
  allocate( &
    leftpressure(Nface), rightpressure(Nface), deltap(Nface), &
    leftparameter(Nface), rightparameter(Nface), &
    leftdepth(Nface), rightdepth(Nface), &
    suliciuvelocity(Nface), suliciupressure(Nface), &
    leftvelocity(Nface), rightvelocity(Nface), &
    fluxleft(2,Nface), fluxright(2,Nface) &
  )
  !do iface = 1,Nface !! to add choice other than Suliciu relaxation
  !end do
    ! suliciu parameter initialization
    leftparameter = sqrt(g*cell0(1:N0-1)%depth)*cell0(1:N0-1)%depth
    rightparameter = sqrt(g*cell0(2:N0)%depth)*cell0(2:N0)%depth
    ! initialization satisfying subcharacteristic condition
    call suliciu_initialization( leftparameter, rightparameter, cell0 ) 
    ! intermediate (suliciu) velocity computation (riemann invariant)
    suliciuvelocity = ( leftparameter*cell0(1:N0-1)%velocity + rightparameter*cell0(2:N0)%velocity &
      + cell0(1:N0-1)%pressure - cell0(2:N0)%pressure )/(rightparameter+leftparameter)
    ! intermediate (suliciu) pressure computation (riemann invariant)
    suliciupressure = ( cell0(2:N0)%pressure*leftparameter + cell0(1:N0-1)%pressure*rightparameter &
      + leftparameter*rightparameter*(cell0(1:N0-1)%velocity-cell0(2:N0)%velocity) )/(rightparameter+leftparameter)
    ! intermediate (suliciu) depths
    leftdepth = 1./( 1./cell0(1:N0-1)%depth + (cell0(1:N0-1)%pressure-suliciupressure)/leftparameter**2 )
    rightdepth = 1./( 1./cell0(2:N0)%depth + (cell0(2:N0)%pressure-suliciupressure)/rightparameter**2 )
    ! intermediate (suliciu) eigenvalue computation (contact discontinuity velocity)
    !leftvelocity = cell0(1:N0-1)%velocity - leftparameter/cell0(1:N0-1)%depth
    !rightvelocity = cell0(2:N0)%velocity + rightparameter/cell0(2:N0)%depth
    leftvelocity = suliciuvelocity - leftparameter/leftdepth
    rightvelocity = suliciuvelocity + rightparameter/rightdepth
    ! flux computation
    where( (suliciuvelocity>=0.).AND.( leftvelocity>=0) )
      fluxleft(1,:) = cell0(1:N0-1)%depth*cell0(1:N0-1)%velocity
      fluxleft(2,:) = cell0(1:N0-1)%depth*cell0(1:N0-1)%velocity**2 + cell0(1:N0-1)%pressure
      fluxright(1,:) = fluxleft(1,:)
      fluxright(2,:) = fluxleft(2,:)
    end where
    where( (suliciuvelocity>=0.).AND.( leftvelocity<0) )
      fluxleft(1,:) = leftdepth*suliciuvelocity
      fluxleft(2,:) = leftdepth*suliciuvelocity**2 + suliciupressure
      fluxright(1,:) = fluxleft(1,:)
      fluxright(2,:) = fluxleft(2,:)
    end where
    where( (suliciuvelocity<0.).AND.( rightvelocity>=0) )
      fluxleft(1,:) = rightdepth*suliciuvelocity
      fluxleft(2,:) = rightdepth*suliciuvelocity**2 + suliciupressure
      fluxright(1,:) = fluxleft(1,:)
      fluxright(2,:) = fluxleft(2,:)
    end where
    where( (suliciuvelocity<0.).AND.( rightvelocity<0) )
      fluxleft(1,:) = cell0(2:N0)%depth*cell0(2:N0)%velocity
      fluxleft(2,:) = cell0(2:N0)%depth*cell0(2:N0)%velocity**2 + cell0(2:N0)%pressure 
      fluxright(1,:) = fluxleft(1,:)
      fluxright(2,:) = fluxleft(2,:)
    end where
  where( rightvelocity<0 ) rightvelocity = -rightvelocity
  where( leftvelocity<0 ) leftvelocity = -leftvelocity
  !dt = CFL*dx/max(maxval(rightvelocity),maxval(leftvelocity))
  dt = CFL*minval( (/cell0(2:N0)%volume/rightvelocity,cell0(1:N0-1)%volume/leftvelocity/) )
  if (dt<dtmin) then
    dt = dtmin
    print *, 'dt too small'
  else
    dt = min(dt,dtmax) 
  end if
  cell0(2:N0-1)%depth = cell0(2:N0-1)%depth + dt*fluxleft(1,1:Nface-1)/cell0(2:N0-1)%volume
  cell0(2:N0-1)%depth = cell0(2:N0-1)%depth - dt*fluxright(1,2:Nface)/cell0(2:N0-1)%volume 
  cell0(2:N0-1)%discharge = cell0(2:N0-1)%discharge + dt*fluxleft(2,1:Nface-1)/cell0(2:N0-1)%volume
  cell0(2:N0-1)%discharge = cell0(2:N0-1)%discharge - dt*fluxright(2,2:Nface)/cell0(2:N0-1)%volume 
  cell0(2:N0-1)%velocity = cell0(2:N0-1)%discharge/cell0(2:N0-1)%depth
  cell0(2:N0-1)%pressure = g*cell0(2:N0-1)%depth**2 
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  subroutine suliciu_initialization( leftparameter, rightparameter, cell0 ) 
    use m_cell
    implicit none
    double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
    type(t_cell), dimension(:), intent(in) :: cell0
    integer :: N0, Nface
    double precision, dimension(:), allocatable :: temp1, temp2, temp3
    double precision :: alpha
    N0 = size(cell0)
    Nface = N0-1
    allocate( temp1(Nface), temp2(Nface), temp3(Nface) )
    alpha = 2. !1.5
    temp1 = cell0(2:N0)%pressure-cell0(1:N0-1)%pressure
    where(temp1>=0.)
      temp2 = temp1/(rightparameter+leftparameter)
      temp3 = 0.
    elsewhere
      temp2 = 0.
      temp3 = -temp1/(rightparameter+leftparameter)
    end where
    temp1 = cell0(1:N0-1)%velocity-cell0(2:N0)%velocity
    where(temp1<0.) temp1=0.
    temp2 = ( temp2 + temp1 )*cell0(1:N0-1)%depth/leftparameter
    temp3 = ( temp3 + temp1 )*cell0(2:N0)%depth/rightparameter
    ! --
    !leftparameter = leftparameter * ( .5 + temp2 + sqrt( .25 + temp2 ) )
    !rightparameter = rightparameter * ( .5 + temp3 + sqrt( .25 + temp3 ) )
    leftparameter = leftparameter * ( 1. + alpha*temp2 )
    rightparameter = rightparameter * ( 1. + alpha*temp3 )
  end subroutine
end subroutine

! ------------------------------------------------------------

subroutine printout(stencil)
  use m_data
  implicit none
  integer :: stencil, ix
  character (len=20) :: mystring ! buffer string
  print *, 'Cell values including boundary ghost cells'
  write(mystring,*) stencil
  print '('//mystring//'f6.2,"   ...",'//mystring//'f6.2)', &
    (/ (cell(ix)%depth, ix=1,stencil) , &
       (cell(ix)%depth, ix=Nx+2-stencil+1,Nx+2) /) 
end subroutine printout

! subroutine writeout(stencil)
!   use m_data
!   implicit none
!   integer :: stencil, ix
!   character (len=20) :: mystring ! buffer string
! end subroutine printout
