program sv
  use m_data
  implicit none
  !--------------------------------------------------------------------------
  integer :: Neq ! number of equations to be solved
  double precision :: dt    ! current time step
  integer :: nt = 0            ! time iteration number
  double precision :: t = 0.   ! current time
  integer :: ix                ! spatial DOF (cells) iterator
  type (t_cell), dimension(:), allocatable :: cell1 ! buffer cell
  double precision, dimension(:,:), allocatable :: fluxleft, fluxright
  double precision :: t_neighbour ! maximal time allowed by CFL=1
  character (len=20) :: mystring ! buffer string
  logical :: mylogical
  integer :: myinteger, stencil
  double precision :: myzeromachine = 1e-16
  !--------------------------------------------------------------------------
  interface
    subroutine riemann(fluxleft,fluxright,t_neighbour,cell0,g)
      use m_cell
      double precision :: g
      type (t_cell), dimension(:) :: cell0 
      double precision :: t_neighbour
      double precision, dimension(:,:) :: fluxleft, fluxright
    end subroutine riemann
    subroutine printout(stencil,cell0)
      use m_cell
      integer :: stencil
      type (t_cell), dimension(:) :: cell0 
    end subroutine printout
    subroutine writeout(cell)
      use m_cell
      type (t_cell), dimension(:) :: cell 
    end subroutine writeout
  end interface
  ! Initialization ----------------------------------------------------------------
  print *, 'Initialization'
  call svini ! <<<<<<<<<<<<<<<<<<<<<<<< INITIALIZE THE (GLOBAL) VARIABLES OF m_data
  print *, 'The step size in space is ', dx
  print *, 'Post-processing every ', dt_clock
  stencil = 5 ! (Nx+2)/2 ! number of cells shown on left and right boundaries <= (Nx+2)/2
  if(viscoelastic==.FALSE.) Neq = 3
  if(viscoelastic==.TRUE.) Neq = 5
  allocate( cell1(Nx+2), fluxleft(Neq,Nx+1), fluxright(Neq,Nx+1) )
  ! Parameterization ---------------------------------------------------------------
  open(unit=0,file='dt.res',form='formatted',status='new')
  open(unit=1,file='h.res',form='formatted',status='new')
  open(unit=2,file='Q.res',form='formatted',status='new')
  open(unit=3,file='u.res',form='formatted',status='new')
  open(unit=4,file='P.res',form='formatted',status='new')
  open(unit=5,file='phi.res',form='formatted',status='new')
  open(unit=6,file='sigmaxx.res',form='formatted',status='new')
  open(unit=7,file='sigmazz.res',form='formatted',status='new')
  ! Post-processing and boundary conditions applied to copy cell1 -------------------
  cell1 = cell
  !cell1 = trunc(cell,myzeromachine) ! to handle vacuum in Riemann solver
  cell1(1) = cell1(2) ! no-flux
  cell1(Nx+2) = cell1(Nx+1) ! no-flux
  call printout(stencil,cell1)
  print '("* Saving data at ",f6.2," >= ",f6.2)', t, t_clock
  call writeout(cell1)
  t_clock = dt_clock ! next post-processing
  ! Computations -----------------------------------------------------------------
  print *, 'Entering time loop'
  time_loop : do while((t<Tmax).AND.(nt<Ntmax))
    nt = nt+1 
    ! Flux obtained from homogeneous (approximate) Riemann problems --------------
    call riemann(fluxleft,fluxright,t_neighbour,cell1,g) !! Suliciu relaxation : choice ?
    cell1(2:Nx+1)%depth = (fluxleft(1,1:Nx)-fluxright(1,2:Nx+1))/cell(2:Nx+1)%volume
    cell1(2:Nx+1)%discharge = (fluxleft(2,1:Nx)-fluxright(2,2:Nx+1))/cell(2:Nx+1)%volume
    cell1(2:Nx+1)%tracer = (fluxleft(3,1:Nx)-fluxright(3,2:Nx+1))/cell(2:Nx+1)%volume
    if(viscoelastic==.TRUE.) then
      cell1(2:Nx+1)%sigmaxx = (fluxleft(4,1:Nx)-fluxright(4,2:Nx+1))/cell(2:Nx+1)%volume
      cell1(2:Nx+1)%sigmazz = (fluxleft(5,1:Nx)-fluxright(5,2:Nx+1))/cell(2:Nx+1)%volume
    end if
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
      !dt = minval( (\ Tmax-t, dt_clock \) ) 
      dt = min(Tmax-t,t_clock-t) 
    end if
    print '("Time step at iteration",i10," = ",e10.3," time = ",f10.2)',nt,dt,t
    thist(1,nt) = t
    thist(2,nt) = t_neighbour*CFL
    thist(3,nt) = dt
    write( 0, '(3e15.6)') thist(1:3,nt)
    t = t+dt
    ! First-order time-splitting ------------------------------------------------
    cell(2:Nx+1)%depth = cell(2:Nx+1)%depth + dt*cell1(2:Nx+1)%depth
    cell(2:Nx+1)%discharge = cell(2:Nx+1)%discharge + dt*cell1(2:Nx+1)%discharge
    if(minval(cell(2:Nx+1)%depth)<0.) then
      print *, 'Negative depth!'
      exit
    end if
    !cell(2:Nx+1)%velocity = cell1(2:Nx+1)%discharge/cell1(2:Nx+1)%depth !! /0
    where( cell(2:Nx+1)%depth>0. )
      cell(2:Nx+1)%velocity = cell(2:Nx+1)%discharge/cell(2:Nx+1)%depth
    elsewhere
      cell(2:Nx+1)%velocity = 0.
    end where
    cell(2:Nx+1)%tracer = cell(2:Nx+1)%tracer + dt*cell1(2:Nx+1)%tracer
    cell(2:Nx+1)%pressure = g*cell(2:Nx+1)%depth**2/2
    if(viscoelastic==.TRUE.) then
      cell(2:Nx+1)%sigmaxx = cell(2:Nx+1)%sigmaxx + dt*cell1(2:Nx+1)%sigmaxx
      cell(2:Nx+1)%sigmazz = cell(2:Nx+1)%sigmazz + dt*cell1(2:Nx+1)%sigmazz
    end if
    ! Post-processing and boundary conditions applied to copy cell1 --------------
    cell1 = cell 
    !cell1 = trunc(cell,myzeromachine) ! to handle vacuum in the approximate Riemann solver
    cell1(1) = cell1(2) ! no-flux
    cell1(Nx+2) = cell1(Nx+1) ! no-flux
    if(mylogical) then
      call printout(stencil,trunc(cell1,myzeromachine))
      print '("* Saving data at ",f6.2," >= ",f6.2)', t, t_clock
      call writeout(cell1)
      t_clock = t_clock + dt_clock
    end if
  end do time_loop 
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  function trunc(cell,zeromachine)
    use m_cell
    implicit none
    type (t_cell), dimension(:) :: cell
    type (t_cell), dimension(:), allocatable :: trunc
    double precision :: zeromachine
    allocate( trunc(size(cell)) )
    trunc = cell
    where( cell%depth<=0. ) ! to handle vacuum in the approximate Riemann solver
      trunc%depth=1e-16
      trunc%discharge=0.
      trunc%velocity=0.
      trunc%pressure=0.
      trunc%tracer=0.
    end where
    return
  end function
end program sv

! ------------------------------------------------------------

subroutine riemann(fluxleft,fluxright,t_neighbour,cell0,g)
  use m_cell
  implicit none
  double precision, intent(in) :: g
  type (t_cell), dimension(:), intent(in) :: cell0 
  double precision, intent(out) :: t_neighbour
  double precision, dimension(:,:), intent(out) :: fluxleft, fluxright
  ! Suliciu relaxation solver: diagonal system with additional variables
  integer :: N0, Nface
  double precision, dimension(:), allocatable :: leftspeed, rightspeed  
  double precision, dimension(:), allocatable :: leftparameter, rightparameter  
  double precision, dimension(:), allocatable :: leftdepth, rightdepth
  double precision, dimension(:), allocatable :: suliciuvelocity, suliciupressure
  double precision, dimension(:), allocatable :: leftvelocity, rightvelocity
  logical, dimension(:), allocatable :: mylogical
  integer :: ii = 0
  !--------------------------------------------------------------------------
  ! interface
  !   subroutine suliciu_initialization( leftparameter, rightparameter, cell0 ) 
  !     use m_cell
  !     double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
  !     type(t_cell), dimension(:), intent(in) :: cell0
  !   end subroutine suliciu_initialization
  ! end interface
  !--------------------------------------------------------------------------
  N0 = size(cell0)
  Nface = N0-1
  allocate( &
    leftspeed(Nface), rightspeed(Nface), &
    leftparameter(Nface), rightparameter(Nface), &
    leftdepth(Nface), rightdepth(Nface), &
    suliciuvelocity(Nface), suliciupressure(Nface), &
    leftvelocity(Nface), rightvelocity(Nface), &
    mylogical(Nface) &
  )
  leftspeed = sqrt(g*cell0(1:N0-1)%depth)
  rightspeed = sqrt(g*cell0(2:N0)%depth)
  call suliciu_initialization( leftspeed, rightspeed, cell0 ) 
  leftvelocity = cell0(1:N0-1)%velocity - leftspeed
  rightvelocity = cell0(2:N0)%velocity + rightspeed
  !!! suliciu Lagrangian speeds initialization satisfying subcharacteristic condition
  leftparameter = leftspeed*cell0(1:N0-1)%depth
  rightparameter = rightspeed*cell0(2:N0)%depth
  !where( (cell0(1:N0-1)%depth==0.).AND.(cell0(2:N0)%depth==0.) )
    suliciuvelocity = 0.
    suliciupressure = 0.
    leftdepth = 0.
    rightdepth = 0.
  !end where
  where( (cell0(1:N0-1)%depth/=0.).AND.(cell0(2:N0)%depth==0.) )
    suliciuvelocity = cell0(1:N0-1)%velocity + cell0(1:N0-1)%pressure/leftparameter
    !suliciupressure = 0.
    leftdepth = cell0(1:N0-1)%depth/( 1. + cell0(1:N0-1)%pressure/(leftspeed*leftparameter) )
    !rightdepth = 0.
  end where
  where( (cell0(1:N0-1)%depth==0.).AND.(cell0(2:N0)%depth/=0.) )
    suliciuvelocity = cell0(2:N0)%velocity - cell0(2:N0)%pressure/rightparameter
    !suliciupressure = 0.
    !leftdepth = 0.
    rightdepth = cell0(2:N0)%depth/( 1. + cell0(2:N0)%pressure/(rightspeed*rightparameter) )
  end where
  where( (cell0(1:N0-1)%depth/=0.).AND.(cell0(2:N0)%depth/=0.) )
    suliciuvelocity = ( leftparameter*cell0(1:N0-1)%velocity + rightparameter*cell0(2:N0)%velocity &
      + cell0(1:N0-1)%pressure - cell0(2:N0)%pressure )/(rightparameter+leftparameter)
    suliciupressure = ( cell0(2:N0)%pressure*leftparameter + cell0(1:N0-1)%pressure*rightparameter &
      + leftparameter*rightparameter*(cell0(1:N0-1)%velocity-cell0(2:N0)%velocity) )/(rightparameter+leftparameter)
    leftdepth = cell0(1:N0-1)%depth / &
    ( 1. - (( suliciupressure - cell0(1:N0-1)%pressure ) / ( leftspeed*leftparameter )) )
    rightdepth = cell0(2:N0)%depth / &
    ( 1. - (( suliciupressure - cell0(2:N0)%pressure ) / ( rightspeed*rightparameter )) )
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
!   mylogical = (cell0(1:N0-1)%depth==0.).AND.(cell0(2:N0)%depth==0.)
!   ii = 1
!   do while( mylogical(ii).eqv..FALSE. )
!     ii = ii+1
!   end do
!   print *, ii
!   print *, cell0(ii-1)%depth
!   print *, cell0(ii)%depth
!   print '(20e10.3)', cell0%depth
!   print '(20f10.3)', cell0%velocity
!   print '(20f10.3)', cell0%pressure
!   print '(20e10.3)', leftspeed
!   print '(20e10.3)', rightspeed
!   print '(21f10.3)', leftvelocity
!   print '(21f10.3)', suliciuvelocity
!   print '(21f10.3)', rightvelocity
!   print '(21f10.3)', suliciupressure
!   print '(20e10.3)', leftdepth
!   print '(20e10.3)', rightdepth
!   print '(21f10.3)', leftparameter
!   print '(21f10.3)', rightparameter
!   print *, '---------DEBUG-------------'
  !!! flux computation
  if(viscoelastic==.FALSE.) then
  where( (suliciuvelocity>=0.).AND.(leftvelocity>=0.) )
    fluxleft(1,:) = cell0(1:N0-1)%depth*cell0(1:N0-1)%velocity
    fluxleft(2,:) = cell0(1:N0-1)%depth*cell0(1:N0-1)%velocity**2 + cell0(1:N0-1)%pressure
    fluxleft(3,:) = cell0(1:N0-1)%tracer*cell0(1:N0-1)%velocity
  end where
  where( (suliciuvelocity>=0.).AND.(leftvelocity<0.) )
    fluxleft(1,:) = leftdepth*suliciuvelocity
    fluxleft(2,:) = leftdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = cell0(1:N0-1)%tracer*suliciuvelocity
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity>=0.) )
    fluxleft(1,:) = rightdepth*suliciuvelocity
    fluxleft(2,:) = rightdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = cell0(2:N0)%tracer*suliciuvelocity
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity<0.) )
    fluxleft(1,:) = cell0(2:N0)%depth*cell0(2:N0)%velocity
    fluxleft(2,:) = cell0(2:N0)%depth*cell0(2:N0)%velocity**2 + cell0(2:N0)%pressure 
    fluxleft(3,:) = cell0(2:N0)%tracer*cell0(2:N0)%velocity
  end where
  fluxright = fluxleft ! conservative
  end if
  if(viscoelastic==.TRUE.) then
  where( (suliciuvelocity>=0.).AND.(leftvelocity>=0.) )
    fluxleft(1,:) = cell0(1:N0-1)%depth*cell0(1:N0-1)%velocity
    fluxleft(2,:) = cell0(1:N0-1)%depth*cell0(1:N0-1)%velocity**2 + cell0(1:N0-1)%pressure
    fluxleft(3,:) = cell0(1:N0-1)%tracer*cell0(1:N0-1)%velocity
    fluxleft(4,:) = cell0(1:N0-1)%sigmaxx*cell0(1:N0-1)%discharge
    fluxleft(5,:) = cell0(1:N0-1)%sigmazz*cell0(1:N0-1)%discharge
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
    fluxleft(4,:) = cell0(2:N0)%sigmaxx*cell0(2:N0)%discharge & 
      - leftvelocity*
      - suliciuvelocity*
      - rightvelocity*
    fluxleft(5,:) = cell0(2:N0)%sigmazz*cell0(2:N0)%discharge
  end where
  where( (suliciuvelocity>=0.).AND.(leftvelocity<0.) )
    fluxleft(1,:) = leftdepth*suliciuvelocity
    fluxleft(2,:) = leftdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = cell0(1:N0-1)%tracer*suliciuvelocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity>=0.) )
    fluxleft(1,:) = rightdepth*suliciuvelocity
    fluxleft(2,:) = rightdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = cell0(2:N0)%tracer*suliciuvelocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity<0.) )
    fluxleft(1,:) = cell0(2:N0)%depth*cell0(2:N0)%velocity
    fluxleft(2,:) = cell0(2:N0)%depth*cell0(2:N0)%velocity**2 + cell0(2:N0)%pressure 
    fluxleft(3,:) = cell0(2:N0)%tracer*cell0(2:N0)%velocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
  end where
  end if
!   print *, '---------DEBUG-------------'
!   print '(21f10.3)', fluxleft(1,:)
!   print '(21f10.3)', fluxleft(2,:)
!   print '(21f10.3)', fluxleft(3,:)
!   print *, '---------DEBUG-------------'
  where( rightvelocity<0 ) rightvelocity = -rightvelocity
  where( leftvelocity<0 ) leftvelocity = -leftvelocity
  t_neighbour = minval( (/cell0(2:N0)%volume/rightvelocity,cell0(1:N0-1)%volume/leftvelocity/) )
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  subroutine suliciu_initialization( leftparameter, rightparameter, cell0 ) 
    use m_cell
    implicit none
    double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
    type(t_cell), dimension(:), intent(in) :: cell0
    integer :: N0, Nface
    double precision, dimension(:), allocatable :: temp0, temp1, temp2, temp3
    double precision :: alpha
    N0 = size(cell0)
    Nface = N0-1
    allocate( temp0(Nface), temp1(Nface), temp2(Nface), temp3(Nface) )
    temp0 = rightparameter*cell0(2:N0)%depth + leftparameter*cell0(1:N0-1)%depth
    temp1 = cell0(2:N0)%pressure-cell0(1:N0-1)%pressure ! Pr-Pl
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
    temp1 = cell0(1:N0-1)%velocity-cell0(2:N0)%velocity ! u_l-u_r
    where(temp1<0.) temp1=0.
    temp2 = ( temp2 + temp1 )
    temp3 = ( temp3 + temp1 )
    ! --
    alpha = 2. !1.5
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
  print '('//mystring//'f6.2,"   ...",'//mystring//'f6.2)', &
    (/ (cell(ix)%depth, ix=1,stencil) , &
       (cell(ix)%depth, ix=Nx+2-stencil+1,Nx+2) /) 
!   print '('//mystring//'f6.2,"   ...",'//mystring//'f6.2)', &
!     (/ (cell(ix)%tracer, ix=1,stencil) , &
!         (cell(ix)%tracer, ix=Nx+2-stencil+1,Nx+2) /) 
!   print '('//mystring//'f6.2,"   ...",'//mystring//'f6.2)', &
!     (/ (cell(ix)%tracer/cell(ix)%depth, ix=1,stencil) , &
!        (cell(ix)%tracer/cell(ix)%depth, ix=Nx+2-stencil+1,Nx+2) /) 
end subroutine printout

subroutine writeout(cell)
  use m_cell
  implicit none
  type (t_cell), dimension(:) :: cell 
  integer :: Nx
  character (len=20) :: mystring ! buffer string
  Nx = size(cell)
  write(mystring,*) Nx
  write( 1, '('//mystring//'f6.2)') cell(2:(Nx+1))%depth
  write( 2, '('//mystring//'f6.2)') cell(2:(Nx+1))%discharge
  write( 3, '('//mystring//'f6.2)') cell(2:(Nx+1))%velocity
  write( 4, '('//mystring//'f6.2)') cell(2:(Nx+1))%pressure
  write( 5, '('//mystring//'f6.2)') cell(2:(Nx+1))%tracer/ &
    (cell1(2:(Nx+1))%depth+myzeromachine*maxval((cell1(2:(Nx+1))%depth)))
end subroutine printout
