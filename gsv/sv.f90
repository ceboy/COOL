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
    subroutine pressure(cell,theta,gavrilyuk,elasticmodulus,oneoverell)
      use m_physics
      use m_cell
      type (t_cell), dimension(:) :: cell 
      double precision :: theta, gavrilyuk, elasticmodulus, oneoverell
    end subroutine pressure
    subroutine riemann(t_neighbour,interf,alpha)
      use m_cell
      double precision :: t_neighbour
      type (t_cell), dimension(:) :: interf 
      double precision :: alpha
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
  ! Initialization : myzeromachine -------------------------------------------------
  myzeromachine = 1.
  do while(1./=1.+myzeromachine)
    myzeromachine = myzeromachine/2.
  end do
  print '("Zero machine = ",e12.6)', myzeromachine
  print *, 'Initialization'
  ! Initialization : svini ---------------------------------------------------------
  call svini ! <<<<<<<<<<<<<<<<<<<<<<<< INITIALIZE THE (GLOBAL) VARIABLES OF m_data
  print *, 'The step size in space is ', dx
  print *, 'Post-processing every ', dt_clock
  ! Parameterization ---------------------------------------------------------------
  open(unit=0,file='dt.res',form='formatted',status='new')
  open(unit=1,file='h.res',form='formatted',status='new')
  open(unit=2,file='Q.res',form='formatted',status='new')
  open(unit=3,file='u.res',form='formatted',status='new')
  open(unit=10,file='phi.res',form='formatted',status='new')
  open(unit=11,file='sxx.res',form='formatted',status='new')
  open(unit=12,file='szz.res',form='formatted',status='new')
  open(unit=13,file='microphi.res',form='formatted',status='new')
  open(unit=14,file='macrophi.res',form='formatted',status='new')
  !open(unit=4,file='P.res',form='formatted',status='new')
  ! >>> ABOVE: to be improved (writeout dependency on unit numbers)
  ! Post-processing and boundary conditions applied to copy interf -------------------
  allocate( interf(Nx+2) )
  interf = cell ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COPY NOT VERY EFFICIENT
  interf(1) = interf(2) ! no-flux
  interf(Nx+2) = interf(Nx+1) ! no-flux
  stencil = 6 ! number of cells shown on left and right boundaries <= (Nx+2)/2
  call printout(stencil,interf)
  print '("* Saving data at ",f16.8," >= ",f16.8)', t, t_clock
  call writeout(interf)
  t_clock = dt_clock ! next post-processing
  ! Computations -----------------------------------------------------------------
  print *, 'Entering time loop'
  time_loop : do while((t<Tmax).AND.(nt<Ntmax))
    nt = nt+1 
    ! Flux obtained from homogeneous (approximate) Riemann problems --------------
    call riemann(t_neighbour,interf,alphaspeed) !! Suliciu : choice ?
    interf = fluxtrunc(interf) ! <<<< zero machine into account for cleaner (conservative) results
    ! >>>> ABOVE: one should better create a new type "flux" with same attributes as cell !
    ! Time-step ------------------------------------------------------------------
    dt = t_neighbour*CFL
    if (dt<dtmin) then
      print *, 'dt too small' !dt = dtmin
      exit
    end if
    if (oneoverlambda*dt>1.) then
      dt = 1./oneoverlambda
    end if
    mylogical = .FALSE.
    if ((t+dt>Tmax).OR.(t+dt>t_clock)) then
      mylogical = .TRUE. 
      dt = min(Tmax-t,t_clock-t) 
    end if
    !! minval( (/ Tmax-t, t_clock-t, 1./max(myzeromachine,oneoverlambda) /) ) 
    print '("Time step at iteration ",i10," = ",e10.3," time = ",f10.2)',nt,dt,t
    thist(1,nt) = t ; thist(2,nt) = t_neighbour*CFL ; thist(3,nt) = dt
    write( 0, '(3e15.6)') thist(1:3,nt)
    t = t+dt
    ! First-order time-splitting: step 1 fluxes --------------------------------
    cell(2:Nx+1)%depth = cell(2:Nx+1)%depth + dt*interf(2:Nx+1)%depth
    if(minval(cell(2:Nx+1)%depth)<0.) then
      print *, 'Negative depth!'; exit
    end if
    cell(2:Nx+1)%discharge = cell(2:Nx+1)%discharge + dt*interf(2:Nx+1)%discharge
    cell(2:Nx+1)%htracer = cell(2:Nx+1)%htracer + dt*interf(2:Nx+1)%htracer
    cell(2:Nx+1)%hsigmaxx = cell(2:Nx+1)%hsigmaxx + dt*interf(2:Nx+1)%hsigmaxx
    cell(2:Nx+1)%hsigmazz = cell(2:Nx+1)%hsigmazz + dt*interf(2:Nx+1)%hsigmazz
    cell(2:Nx+1)%hmicroenstrophy = cell(2:Nx+1)%hmicroenstrophy + dt*interf(2:Nx+1)%hmicroenstrophy
    cell(2:Nx+1)%hmacroenstrophy = cell(2:Nx+1)%hmacroenstrophy + dt*interf(2:Nx+1)%hmacroenstrophy
    ! First-order time-splitting: step 2 dissipative sources --------------------
    cell(2:Nx+1)%hsigmaxx = &
      (1.-dt*oneoverlambda)*cell(2:Nx+1)%hsigmaxx + dt*oneoverlambda*cell(2:Nx+1)%depth
    cell(2:Nx+1)%hsigmazz = &
      (1.-dt*oneoverlambda)*cell(2:Nx+1)%hsigmazz + dt*oneoverlambda*cell(2:Nx+1)%depth
    ! First-order time-splitting: step 3 forcing sources --------------------
    cell(2:Nx+1)%discharge = cell(2:Nx+1)%discharge + dt*(g*tan(theta))*cell(2:Nx+1)%depth
    ! Post-processing: output + prepare next step -------------------------------
    cell = celltrunc(cell) ! <<<< zero machine into account for cleaner (extensive) results
    where( (cell(2:Nx+1)%depth>0.) )
      cell(2:Nx+1)%velocity = cell(2:Nx+1)%discharge/cell(2:Nx+1)%depth
      cell(2:Nx+1)%tracer = cell(2:Nx+1)%htracer/cell(2:Nx+1)%depth
      cell(2:Nx+1)%sigmaxx = cell(2:Nx+1)%hsigmaxx/cell(2:Nx+1)%depth
      cell(2:Nx+1)%sigmazz = cell(2:Nx+1)%hsigmazz/cell(2:Nx+1)%depth
      cell(2:Nx+1)%microenstrophy = cell(2:Nx+1)%hmicroenstrophy/cell(2:Nx+1)%depth
      cell(2:Nx+1)%macroenstrophy = cell(2:Nx+1)%hmacroenstrophy/cell(2:Nx+1)%depth
    elsewhere
      cell(2:Nx+1)%velocity = 0.
      cell(2:Nx+1)%tracer = 0.
      cell(2:Nx+1)%sigmaxx = 0.
      cell(2:Nx+1)%sigmazz = 0.
      cell(2:Nx+1)%microenstrophy = 0.
      cell(2:Nx+1)%macroenstrophy = 0.
    end where
    ! First-order time-splitting: step 2 dissipative sources --------------------
!     cell(2:Nx+1)%hsigmaxx = &
!       ((1.-dt*oneoverlambda)*cell(2:Nx+1)%sigmaxx + dt*oneoverlambda)*cell(2:Nx+1)%depth
!     cell(2:Nx+1)%hsigmazz = &
!       ((1.-dt*oneoverlambda)*cell(2:Nx+1)%sigmazz + dt*oneoverlambda)*cell(2:Nx+1)%depth
    !
    call pressure(cell,theta,gavrilyuk,elasticmodulus,oneoverell)
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
      celltrunc%hmicroenstrophy = 0.
      celltrunc%hmacroenstrophy = 0.
      celltrunc%velocity = 0.
      celltrunc%tracer = 0.
      celltrunc%sigmaxx = 0.
      celltrunc%sigmazz = 0.
      celltrunc%microenstrophy = 0.
      celltrunc%macroenstrophy = 0.
      celltrunc%pressure = 0.
      celltrunc%speed = 0.
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
    !
    ref = max(maxval(interf%depth),-minval(interf%depth))
    where( (interf%depth<=myzeromachine*ref).AND.(interf%depth>0.) ) fluxtrunc%depth=0.
    where( (interf%depth>=-myzeromachine*ref).AND.(interf%depth<0.) ) fluxtrunc%depth=0.
    !
    ref = max(maxval(interf%discharge),-minval(interf%discharge))
    where( (interf%discharge<=myzeromachine*ref).AND.(interf%discharge>0.) ) fluxtrunc%discharge=0.
    where( (interf%discharge>=-myzeromachine*ref).AND.(interf%discharge<0.) ) fluxtrunc%discharge=0.
    !
    ref = max(maxval(interf%htracer),-minval(interf%htracer))
    where( (interf%htracer<=myzeromachine*ref).AND.(interf%htracer>0.) ) fluxtrunc%htracer=0.
    where( (interf%htracer>=-myzeromachine*ref).AND.(interf%htracer<0.) ) fluxtrunc%htracer=0.
    !
    ref = max(maxval(interf%hsigmaxx),-minval(interf%hsigmaxx))
    where( (interf%hsigmaxx<=myzeromachine*ref).AND.(interf%hsigmaxx>0.) ) fluxtrunc%hsigmaxx=0.
    where( (interf%hsigmaxx>=-myzeromachine*ref).AND.(interf%hsigmaxx<0.) ) fluxtrunc%hsigmaxx=0.
    !
    ref = max(maxval(interf%hsigmazz),-minval(interf%hsigmazz))
    where( (interf%hsigmazz<=myzeromachine*ref).AND.(interf%hsigmazz>0.) ) fluxtrunc%hsigmazz=0.
    where( (interf%hsigmazz>=-myzeromachine*ref).AND.(interf%hsigmazz<0.) ) fluxtrunc%hsigmazz=0.
    !
    ref = max(maxval(interf%hmicroenstrophy),-minval(interf%hmicroenstrophy))
    where( (interf%hmicroenstrophy<=myzeromachine*ref).AND.(interf%hmicroenstrophy>0.) )     fluxtrunc%hmicroenstrophy=0.
    where( (interf%hmicroenstrophy>=-myzeromachine*ref).AND.(interf%hmicroenstrophy<0.) ) fluxtrunc%hmicroenstrophy=0.
    !
    ref = max(maxval(interf%hmacroenstrophy),-minval(interf%hmacroenstrophy))
    where( (interf%hmacroenstrophy<=myzeromachine*ref).AND.(interf%hmacroenstrophy>0.) ) fluxtrunc%hmacroenstrophy=0.
    where( (interf%hmacroenstrophy>=-myzeromachine*ref).AND.(interf%hmacroenstrophy<0.) ) fluxtrunc%hmacroenstrophy=0.
    return
  end function
end program sv

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
  print '('//mystring//'f16.8,"   ...",'//mystring//'f16.8)', &
    (/ (cell(ix)%pressure, ix=1,stencil) , &
        (cell(ix)%pressure, ix=Nx+2-stencil+1,Nx+2) /) 
!   print '('//mystring//'f16.8,"   ...",'//mystring//'f16.8)', &
!     (/ (cell(ix)%sigmazz, ix=1,stencil) , &
!        (cell(ix)%sigmazz, ix=Nx+2-stencil+1,Nx+2) /) 
end subroutine printout

! ------------------------------------------------------------

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
  write( 10, '('//mystring//'f16.8)') cell(2:(Nx+1))%tracer
  write( 11, '('//mystring//'f16.8)') cell(2:(Nx+1))%sigmaxx
  write( 12, '('//mystring//'f16.8)') cell(2:(Nx+1))%sigmazz
  write( 13, '('//mystring//'f16.8)') cell(2:(Nx+1))%microenstrophy
  write( 14, '('//mystring//'f16.8)') cell(2:(Nx+1))%macroenstrophy
  !write( 4, '('//mystring//'f16.8)') cell(2:(Nx+1))%pressure
end subroutine writeout


