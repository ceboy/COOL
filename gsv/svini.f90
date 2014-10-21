subroutine svini
  use m_data
  implicit none
  !--------------------------------------------------------------------------
  integer :: myiarg            ! number of command arguments
  character(len=10) :: Nxstring
  integer :: ix                ! spatial DOF (cells) iterator
  select case(testcase)
  !--------------------------------------------------------------------------
  case(1)
    Tmax = .5d0
    dtmin = 1.d-4
    Ntmax = 100000
    allocate( thist(3,Ntmax) )
    Ntmax_clock = 501 !501
    dt_clock = Tmax/(Ntmax_clock-1)
    !dt_clock = 1.d-1 ! 
    !Ntmax_clock = (Tmax/dt_clock)+1 
    CFL = .5d0
    Nx = 10
    myiarg = iargc() ! compiler dependent ? to read arguments
    print *, myiarg
    if (myiarg>0) then
      call getarg(1,Nxstring)
      read(Nxstring,*) Nx
    end if
    dx = 10.d0/Nx
    allocate( cell(Nx+2) )
    volume_initialization : do ix = 1,Nx+2
      cell(ix)%volume = dx
    end do volume_initialization
    cell(1)%center = -5. - dx/2.
    center_initialization : do ix = 2,Nx+1
      cell(ix)%center = cell(ix-1)%center + (cell(ix-1)%volume+cell(ix)%volume)/2.
    end do center_initialization
    depth_initialization : do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%depth = 3.d0 
      else
	cell(ix)%depth = 1.d0
      endif
    end do depth_initialization
    discharge_initialization : do ix = 1,Nx+2
      cell(ix)%discharge = 0d0
    end do discharge_initialization
    velocity_initialization : do ix = 1,Nx+2
      cell(ix)%velocity = 0d0
    end do velocity_initialization
    pressure_initialization : do ix = 1,Nx+2
      cell(ix)%pressure = g*cell(ix)%depth**2
    end do pressure_initialization
    tracer_initialization : do ix = 1,Nx+2
      if (cell(ix)%center<0.) then
	cell(ix)%tracer = 2.d0 
      else
	cell(ix)%tracer = 1.d0
      endif
    end do tracer_initialization
  !--------------------------------------------------------------------------
  end select
end subroutine svini