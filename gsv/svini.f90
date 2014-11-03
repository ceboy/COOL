subroutine svini
  use m_physics
  use m_data
  implicit none
  !--------------------------------------------------------------------------
  integer :: myiarg            ! number of command arguments
  character(len=10) :: Nxstring
  integer :: ix                ! spatial DOF (cells) iterator
  select case(testcase)
  !--------------------------------------------------------------------------
  case(1)
    elasticmodulus = 0
    Tmax = .1d0
    dtmin = 1.d-6
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%depth = 3.d0 
      else
	cell(ix)%depth = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0.d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0.d0
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%tracer = 2.d0
      else
	cell(ix)%tracer = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2.
      cell(ix)%speed = sqrt( g*cell(ix)%depth )
    end do 
  !--------------------------------------------------------------------------
  case(2)
    elasticmodulus = 0
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%depth = 3.d0 
      else
	cell(ix)%depth = 0.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%tracer = 2.d0
      else
	cell(ix)%tracer = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2.
      cell(ix)%speed = sqrt( g*cell(ix)%depth )
    end do 
  !--------------------------------------------------------------------------
  case(3)
    elasticmodulus = 0
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center<1.).AND.(cell(ix)%center>-1.)) then
	cell(ix)%depth = 3.d0 
      else
	cell(ix)%depth = 0.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center<.5).AND.(cell(ix)%center>-.5)) then
	cell(ix)%tracer = 2.d0
      else
	cell(ix)%tracer = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2.
      cell(ix)%speed = sqrt( g*cell(ix)%depth )
    end do 
  !--------------------------------------------------------------------------
  case(4)
    elasticmodulus = 0
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center>1.).OR.(cell(ix)%center<-1.)) then
	cell(ix)%depth = 5.d0 
      else
	cell(ix)%depth = 0.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center>2.).OR.(cell(ix)%center<-2.)) then
	cell(ix)%tracer = 2.d0
      else
	cell(ix)%tracer = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2.
      cell(ix)%speed = sqrt( g*cell(ix)%depth )
    end do 
  !--------------------------------------------------------------------------
  case(5)
    elasticmodulus = 0
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center<1.).AND.(cell(ix)%center>-1.)) then
	cell(ix)%depth = 2.d0 
      else if (cell(ix)%center>1.) then
	cell(ix)%depth = 1.d0
      else if (cell(ix)%center<-1.) then
	cell(ix)%depth = 10.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center>-1.).AND.(cell(ix)%center<1.)) then
	cell(ix)%tracer = 1.d0
      else
	cell(ix)%tracer = 0.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2.
      cell(ix)%speed = sqrt( g*cell(ix)%depth )
    end do 
  !--------------------------------------------------------------------------
  case(6)
    elasticmodulus = 0
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center<1.).AND.(cell(ix)%center>-1.)) then
	cell(ix)%depth = 1.d0 + sin(kwave*Pi*(cell(ix)%center+1.))
      else if (cell(ix)%center>1.) then
	cell(ix)%depth = 1.d0
      else if (cell(ix)%center<-1.) then
	cell(ix)%depth = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if ((cell(ix)%center>-1.).AND.(cell(ix)%center<1.)) then
	cell(ix)%tracer = 1.d0
      else
	cell(ix)%tracer = 0.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2.
      cell(ix)%speed = sqrt( g*cell(ix)%depth )
    end do 
  !--------------------------------------------------------------------------
  case(7)
    elasticmodulus = 0
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%depth = 3.d0 
      else
	cell(ix)%depth = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%tracer = 2.d0
      else
	cell(ix)%tracer = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%sigmaxx = 1.d0
    end do
    do ix = 1,Nx+2
      cell(ix)%sigmazz = 1.d0
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2. & 
        + elasticmodulus*(cell(ix)%sigmazz-cell(ix)%sigmaxx)
      cell(ix)%speed = sqrt( g*cell(ix)%depth &
        + elasticmodulus*(3*cell(ix)%sigmazz+cell(ix)%sigmaxx) )
    end do 
  !--------------------------------------------------------------------------
  case(8)
    elasticmodulus = 10.
    Tmax = .5d0
    dtmin = 1.d-6
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
        cell(ix)%depth = 3.d0 
      else
        cell(ix)%depth = .000d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
        cell(ix)%tracer = 2.d0
        cell(ix)%sigmaxx = 1.d0
        cell(ix)%sigmazz = 1.d0
      else
        cell(ix)%tracer = 0.d0
        cell(ix)%sigmaxx = 0.d0
        cell(ix)%sigmazz = 0.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2. & 
        + elasticmodulus*(cell(ix)%sigmazz-cell(ix)%sigmaxx)
      cell(ix)%speed = sqrt( g*cell(ix)%depth &
        + elasticmodulus*(3*cell(ix)%sigmazz+cell(ix)%sigmaxx) )
    end do 
  !--------------------------------------------------------------------------
  case(9)
    elasticmodulus = 0.
    oneoverell = 0. !1./4.
    Tmax = .5d0
    dtmin = 1.d-6
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 51 ; dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-2 ; Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do
    cell(1)%center = -5. - dx/2.
    do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
        cell(ix)%depth = 3.d0 
      else
        cell(ix)%depth = 0.000000d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do
    do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do
    do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
        cell(ix)%tracer = 2.d0
        cell(ix)%sigmaxx = 1.d0
        cell(ix)%sigmazz = 1.d0
      else
        cell(ix)%tracer = 1.d0
        cell(ix)%sigmaxx = 1.d0
        cell(ix)%sigmazz = 1.d0
      endif
    end do
    do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2/2 &
	+ elasticmodulus*(cell(ix)%hsigmazz-cell(ix)%hsigmaxx)/ &
	    (1 + oneoverell*(cell(ix)%sigmazz+cell(ix)%sigmaxx))
      cell(ix)%speed = sqrt( g*cell(ix)%depth &
	+ elasticmodulus*(3*cell(ix)%sigmazz+cell(ix)%sigmaxx)/ &
	    (1 - oneoverell*(cell(ix)%sigmazz+cell(ix)%sigmaxx)) &
	+ oneoverell*2*elasticmodulus*((cell(ix)%sigmazz-cell(ix)%sigmaxx)/ &
	    (1 - oneoverell*(cell(ix)%sigmazz+cell(ix)%sigmaxx)))**2 )
    end do 
  !--------------------------------------------------------------------------
  end select
  cell%htracer = cell%tracer*cell%depth
  cell%hsigmaxx = cell%sigmaxx*cell%depth
  cell%hsigmazz = cell%sigmazz*cell%depth
end subroutine svini