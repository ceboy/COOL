subroutine riemann(t_neighbour,cell,alphaspeed)
  use m_zeromachine
  use m_cell
  implicit none
  double precision, intent(out) :: t_neighbour
  type (t_cell), dimension(:), intent(inout) :: cell 
  double precision, intent(in) :: alphaspeed
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
  !   subroutine suliciu_initialization( leftparameter, rightparameter, cell, alphaspeed ) 
  !     use m_cell
  !     double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
  !     type(t_cell), dimension(:), intent(in) :: cell
  !     double precision, intent(in) :: alphaspeed
  !   end subroutine suliciu_initialization
  ! end interface
  !--------------------------------------------------------------------------
  N0 = size(cell)
  Nface = N0-1
  allocate( fluxleft(7,Nface), fluxright(7,Nface), & !! << easier way to adjust 7 to attributes ??
    leftspeed(Nface), rightspeed(Nface), &
    leftparameter(Nface), rightparameter(Nface), &
    leftdepth(Nface), rightdepth(Nface), &
    suliciuvelocity(Nface), suliciupressure(Nface), &
    leftvelocity(Nface), rightvelocity(Nface), &
    mylogical(Nface) &
  )
  !!! * explicit initialization of the parameter
  leftspeed = cell(1:N0-1)%speed
  rightspeed = cell(2:N0)%speed
  call suliciu_initialization( leftspeed, rightspeed, cell, alphaspeed )
  !!! * extremal wave speeds of Suliciu system
  leftvelocity = cell(1:N0-1)%velocity - leftspeed ! ==suliciuvelocity-leftspeed*(cell(1:N0-1)%depth/leftdepth)
  rightvelocity = cell(2:N0)%velocity + rightspeed ! ==suliciuvelocity+rightspeed*(cell(2:N0)%depth/rightdepth)
  !!! * intermediate variables
  leftparameter = leftspeed*cell(1:N0-1)%depth
  rightparameter = rightspeed*cell(2:N0)%depth
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
  print *, '---------DEBUG-------------'
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
!   print '(10f10.3)', leftvelocity
!   print '(10f10.3)', cell%velocity
!   print '(10f10.3)', suliciuvelocity
!   print '(10f10.3)', rightvelocity
!   print '(10f10.3)', cell%pressure
!   print '(10f10.3)', suliciupressure
!   print *, cell%depth
!   print *, leftdepth
!   print *, rightdepth
!   print *, cell%hsigmaxx*cell%depth
!   print *, lefthsigmaxx*leftdepth
!   print *, righthsigmaxx*rightdepth
!   print *, cell%hsigmazz/(myzeromachine*maxval(cell%depth)+cell%depth)**3
!   print *, lefthsigmazz/(myzeromachine*maxval(leftdepth)+leftdepth)**3
!   print *, righthsigmazz/(myzeromachine*maxval(rightdepth)+rightdepth)**3
!   print *, minval(cell%hsigmaxx*cell%depth), minval(lefthsigmaxx*leftdepth), minval(righthsigmaxx*rightdepth)
  print '("Max principle <= ",3e16.6)', &
     maxval(cell%hsigmaxx*cell%depth), &
     maxval(lefthsigmaxx*leftdepth), &
     maxval(righthsigmaxx*rightdepth)
! <= max principle
  print '("Min principle >= ",3e16.6)', &   
     minval(cell%hsigmazz/(myzeromachine*maxval(cell%depth)+cell%depth)**3), &
     minval(lefthsigmazz/(myzeromachine*maxval(leftdepth)+leftdepth)**3), &
     minval(righthsigmazz/(myzeromachine*maxval(rightdepth)+rightdepth)**3)
! <= min principle
!   print *, maxval(cell%hsigmazz/(myzeromachine*maxval(cell%depth)+cell%depth)**3), &
!     maxval(lefthsigmazz/(myzeromachine*maxval(leftdepth)+leftdepth)**3), &
!     maxval(righthsigmazz/(myzeromachine*maxval(rightdepth)+rightdepth)**3)
  print *, '---------DEBUG-------------'
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
    fluxleft(6,:) = cell(1:N0-1)%hmicroenstrophy*cell(1:N0-1)%velocity
    fluxleft(7,:) = cell(1:N0-1)%hmacroenstrophy*cell(1:N0-1)%velocity
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
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
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
    fluxleft(6,:) = cell(1:N0-1)%microenstrophy*leftdepth*suliciuvelocity
    fluxleft(7,:) = cell(1:N0-1)%macroenstrophy*leftdepth*suliciuvelocity
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
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
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
    fluxleft(6,:) = cell(2:N0)%microenstrophy*rightdepth*suliciuvelocity
    fluxleft(7,:) = cell(2:N0)%macroenstrophy*rightdepth*suliciuvelocity
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
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
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
    fluxleft(6,:) = cell(2:N0)%hmicroenstrophy*cell(2:N0)%velocity
    fluxleft(7,:) = cell(2:N0)%hmacroenstrophy*cell(2:N0)%velocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = cell(2:N0)%sigmaxx*cell(2:N0)%discharge
!     fluxright(5,:) = cell(2:N0)%sigmazz*cell(2:N0)%discharge
    fluxright(4,:) = cell(2:N0)%hsigmaxx*cell(2:N0)%velocity
    fluxright(5,:) = cell(2:N0)%hsigmazz*cell(2:N0)%velocity
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
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
  cell(2:N0-1)%hmicroenstrophy = (fluxright(6,1:Nface-1)-fluxleft(6,2:Nface))/cell(2:N0-1)%volume
  cell(2:N0-1)%hmacroenstrophy = (fluxright(7,1:Nface-1)-fluxleft(7,2:Nface))/cell(2:N0-1)%volume
  where( rightvelocity<0 ) rightvelocity = -rightvelocity
  where( leftvelocity<0 ) leftvelocity = -leftvelocity
  t_neighbour = minval( (/cell(2:N0)%volume/rightvelocity,cell(1:N0-1)%volume/leftvelocity/) )
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  subroutine suliciu_initialization( leftparameter, rightparameter, cell, alphaspeed ) 
    use m_cell
    implicit none
    double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
    type(t_cell), dimension(:), intent(in) :: cell
    double precision, intent(in) :: alphaspeed
    integer :: N0, Nface
    double precision, dimension(:), allocatable :: temp0, temp1, temp2, temp3
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
    !leftparameter = leftparameter*(1. + alphaspeed*temp2/leftparameter)
    !rightparameter = rightparameter*(1. + alphaspeed*temp3/rightparameter)
    leftparameter = leftparameter + alphaspeed*temp2
    rightparameter = rightparameter + alphaspeed*temp3
  end subroutine suliciu_initialization
end subroutine riemann
