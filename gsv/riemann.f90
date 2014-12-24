subroutine riemann(t_neighbour,leftcells,rightcells,alphaspeed)
  use m_zeromachine
  use m_cell
  implicit none
  double precision, intent(out) :: t_neighbour
  type (t_cell), dimension(:), intent(inout) :: leftcells, rightcells
  double precision, intent(in) :: alphaspeed
  integer :: Nface
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
  Nface = size(leftcells) ! = Ninterf en 1D
  if (size(rightcells)/=Nface) then
    print *, "Error: left and right initial cells do not have the same size in riemann.f90"
    stop
  end if
  allocate( & !! << easier way to adjust 7 to attributes: parameterized derived types ??
    fluxleft(7,Nface), fluxright(7,Nface), &
    leftspeed(Nface), rightspeed(Nface), &
    leftparameter(Nface), rightparameter(Nface), &
    leftdepth(Nface), rightdepth(Nface), &
    suliciuvelocity(Nface), suliciupressure(Nface), &
    leftvelocity(Nface), rightvelocity(Nface), &
    mylogical(Nface) &
  )
  !!! * explicit initialization of the parameter
  leftspeed = leftcells%speed
  rightspeed = rightcells%speed
  call suliciu_initialization( leftspeed, rightspeed, leftcells, rightcells, alphaspeed )
  !!! * extremal wave speeds of Suliciu system
  leftvelocity = leftcells%velocity - leftspeed ! ==suliciuvelocity-leftspeed*(leftcells%depth/leftdepth)
  rightvelocity = rightcells%velocity + rightspeed ! ==suliciuvelocity+rightspeed*(rightcells%depth/rightdepth)
  !!! * intermediate variables
  leftparameter = leftspeed*leftcells%depth
  rightparameter = rightspeed*rightcells%depth
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
  where( (leftcells%depth/=0.).AND.(rightcells%depth==0.) ) 
    suliciuvelocity = leftcells%velocity + leftcells%pressure/leftparameter
    ! suliciupressure = 0.
    leftdepth = leftcells%depth/( 1. + leftcells%pressure/(leftspeed*leftparameter) )
    ! rightdepth = 0.
    lefthsigmaxx = leftcells%hsigmaxx*(1. + leftcells%pressure/(leftspeed*leftparameter))
    lefthsigmazz = leftcells%hsigmazz/(1. + leftcells%pressure/(leftspeed*leftparameter))**3
  end where
  where( (leftcells%depth==0.).AND.(rightcells%depth/=0.) )
    suliciuvelocity = rightcells%velocity - rightcells%pressure/rightparameter
    ! suliciupressure = 0.
    ! leftdepth = 0.
    rightdepth = rightcells%depth/( 1. + rightcells%pressure/(rightspeed*rightparameter) )
    righthsigmaxx = rightcells%hsigmaxx*(1. + rightcells%pressure/(rightspeed*rightparameter))
    righthsigmazz = rightcells%hsigmazz/(1. + rightcells%pressure/(rightspeed*rightparameter))**3
  end where
  where( (leftcells%depth/=0.).AND.(rightcells%depth/=0.) )
!     suliciuvelocity = ( leftspeed*leftcells%discharge + rightspeed*rightcells%discharge &
!                         + leftcells%pressure - rightcells%pressure ) &
!                      /( leftspeed*leftcells%depth + rightspeed*rightcells%depth )
    suliciuvelocity = ( leftspeed*leftcells%discharge + rightspeed*rightcells%discharge &
                        + leftcells%pressure - rightcells%pressure ) &
                     /( rightparameter+leftparameter )
!     suliciuvelocity = ( leftparameter*leftcells%velocity + rightparameter*rightcells%velocity &
!       + leftcells%pressure - rightcells%pressure )/(rightparameter+leftparameter)
!! ABOVE: should we use velocity (e.g. to take into account some truncation ?) ?????
!     suliciupressure = ( rightcells%pressure*leftspeed*leftcells%depth &
!                         + leftcells%pressure*rightspeed*rightcells%depth &
!     + leftspeed*rightspeed*(leftcells%discharge*rightcells%depth-rightcells%discharge*leftcells%depth) ) &
!                      /( leftspeed*leftcells%depth + rightspeed*rightcells%depth )
    suliciupressure = ( rightcells%pressure*leftparameter + leftcells%pressure*rightparameter &
     + leftspeed*rightspeed*(leftcells%discharge*rightcells%depth-rightcells%discharge*leftcells%depth) ) &
                     /( rightparameter+leftparameter )
!     suliciupressure = ( rightcells%pressure*leftparameter + leftcells%pressure*rightparameter &
!       + leftparameter*rightparameter*(leftcells%velocity-rightcells%velocity) )/(rightparameter+leftparameter)
!! ABOVE: should we use velocity (e.g. to take into account some truncation ?) ?????
    leftdepth = leftcells%depth/(1.-((suliciupressure-leftcells%pressure)/(leftspeed*leftparameter)))
    rightdepth = rightcells%depth/(1.-((suliciupressure-rightcells%pressure)/(rightspeed*rightparameter)))
!     leftsigmaxx = leftcells%sigmaxx*(leftcells%depth/leftdepth)**2
!     leftsigmazz = leftcells%sigmazz*(leftdepth/leftcells%depth)**2
!     lefthsigmaxx = leftcells%hsigmaxx*(leftcells%depth/leftdepth)
!     lefthsigmazz = leftcells%hsigmazz*(leftdepth/leftcells%depth)**3
    lefthsigmaxx = leftcells%hsigmaxx*(1.-((suliciupressure-leftcells%pressure)/(leftspeed*leftparameter)))
    lefthsigmazz = leftcells%hsigmazz/(1.-((suliciupressure-leftcells%pressure)/(leftspeed*leftparameter)))**3
!     rightsigmaxx = rightcells%sigmaxx*(rightcells%depth/rightdepth)**2
!     rightsigmazz = rightcells%sigmazz*(rightdepth/rightcells%depth)**2
!     righthsigmaxx = rightcells%hsigmaxx*(rightcells%depth/rightdepth)
!     righthsigmazz = rightcells%hsigmazz*(rightdepth/rightcells%depth)**3
    righthsigmaxx = rightcells%hsigmaxx*(1.-((suliciupressure-rightcells%pressure)/(rightspeed*rightparameter)))
    righthsigmazz = rightcells%hsigmazz/(1.-((suliciupressure-rightcells%pressure)/(rightspeed*rightparameter)))**3
  end where
  print *, '---------DEBUG-------------'
!   mylogical = .FALSE.
!   where( suliciuvelocity>rightvelocity )
!     mylogical = .TRUE.
!   end where
!   where( suliciuvelocity<leftvelocity )
!     mylogical = .TRUE.
!   end where
!   mylogical = (leftcells%depth==0.).AND.(rightcells%depth==0.)
!   ii = 1
!   do while( mylogical(ii).eqv..FALSE. )
!     ii = ii+1
!   end do
!   print *, leftspeed
!   print *, rightspeed
!   print *, leftparameter
!   print *, rightparameter
!   print '(10f10.3)', leftvelocity
!   print '(10f10.3)', suliciuvelocity
!   print '(10f10.3)', rightvelocity
!   print '(10f10.3)', suliciupressure
!   print *, leftdepth
!   print *, rightdepth
!   print *, lefthsigmaxx*leftdepth
!   print *, righthsigmaxx*rightdepth
!   print *, lefthsigmazz/(myzeromachine*maxval(leftdepth)+leftdepth)**3
!   print *, righthsigmazz/(myzeromachine*maxval(rightdepth)+rightdepth)**3
!   print *, minval(lefthsigmaxx*leftdepth), minval(righthsigmaxx*rightdepth)
  print '("Max principle <= ",3e16.6)', &
     maxval(lefthsigmaxx*leftdepth), &
     maxval(righthsigmaxx*rightdepth)
! <= max principle
  print '("Min principle >= ",3e16.6)', &   
     minval(lefthsigmazz/(myzeromachine*maxval(leftdepth)+leftdepth)**3), &
     minval(righthsigmazz/(myzeromachine*maxval(rightdepth)+rightdepth)**3)
! <= min principle
!   print *, maxval(lefthsigmazz/(myzeromachine*maxval(leftdepth)+leftdepth)**3), &
!     maxval(righthsigmazz/(myzeromachine*maxval(rightdepth)+rightdepth)**3)
  print *, '---------DEBUG-------------'
  !!! * Suliciu flux computation
  where( (suliciuvelocity>=0.).AND.(leftvelocity>=0.) )
    fluxleft(1,:) = leftcells%depth*leftcells%velocity
    fluxleft(2,:) = leftcells%depth*leftcells%velocity**2 + leftcells%pressure
!     fluxleft(3,:) = leftcells%tracer*leftcells%discharge
!     fluxleft(4,:) = leftcells%sigmaxx*leftcells%discharge
!     fluxleft(5,:) = leftcells%sigmazz*leftcells%discharge
    fluxleft(3,:) = leftcells%htracer*leftcells%velocity
    fluxleft(4,:) = leftcells%hsigmaxx*leftcells%velocity
    fluxleft(5,:) = leftcells%hsigmazz*leftcells%velocity
    fluxleft(6,:) = leftcells%hmicroenstrophy*leftcells%velocity
    fluxleft(7,:) = leftcells%hmacroenstrophy*leftcells%velocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = rightcells%sigmaxx*rightcells%discharge & 
!       - leftvelocity*(leftsigmaxx*leftdepth-leftcells%sigmaxx*leftcells%depth) &
!       - suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx*leftdepth) &
!       - rightvelocity*(rightcells%sigmaxx*rightcells%depth-rightsigmaxx*rightdepth)
!     fluxright(5,:) = rightcells%sigmazz*rightcells%discharge & 
!       - leftvelocity*(leftsigmazz*leftdepth-leftcells%sigmazz*leftcells%depth) &
!       - suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth) &
!       - rightvelocity*(rightcells%sigmazz*rightcells%depth-rightsigmazz*rightdepth)
    fluxright(4,:) = rightcells%hsigmaxx*rightcells%velocity & 
      - rightvelocity*(rightcells%hsigmaxx-righthsigmaxx) &
      - suliciuvelocity*(righthsigmaxx-lefthsigmaxx) &
      - leftvelocity*(lefthsigmaxx-leftcells%hsigmaxx)
    fluxright(5,:) = rightcells%hsigmazz*rightcells%velocity & 
      - rightvelocity*(rightcells%hsigmazz-righthsigmazz) &
      - suliciuvelocity*(righthsigmazz-lefthsigmazz) &
      - leftvelocity*(lefthsigmazz-leftcells%hsigmazz)
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
  end where
  where( (suliciuvelocity>=0.).AND.(leftvelocity<0.) )
    fluxleft(1,:) = leftdepth*suliciuvelocity
    fluxleft(2,:) = leftdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = leftcells%tracer*leftdepth*suliciuvelocity
!     fluxleft(4,:) = leftcells%sigmaxx*leftcells%discharge &
!       + leftvelocity*(leftsigmaxx*leftdepth-leftcells%sigmaxx*leftcells%depth)
!     fluxleft(5,:) = leftcells%sigmazz*leftcells%discharge &
!       + leftvelocity*(leftsigmazz*leftdepth-leftcells%sigmazz*leftcells%depth)
    fluxleft(4,:) = leftcells%hsigmaxx*leftcells%velocity &
      + leftvelocity*(lefthsigmaxx-leftcells%hsigmaxx)
    fluxleft(5,:) = leftcells%hsigmazz*leftcells%velocity &
      + leftvelocity*(lefthsigmazz-leftcells%hsigmazz)
    fluxleft(6,:) = leftcells%microenstrophy*leftdepth*suliciuvelocity
    fluxleft(7,:) = leftcells%macroenstrophy*leftdepth*suliciuvelocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = rightcells%sigmaxx*rightcells%discharge & 
!       - suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx*leftdepth) &
!       - rightvelocity*(rightcells%sigmaxx*rightcells%depth-rightsigmaxx*rightdepth)
!     fluxright(5,:) = rightcells%sigmazz*rightcells%discharge & 
!       - suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth) &
!       - rightvelocity*(rightcells%sigmazz*rightcells%depth-rightsigmazz*rightdepth)
    fluxright(4,:) = rightcells%hsigmaxx*rightcells%velocity & 
      - rightvelocity*(rightcells%hsigmaxx-righthsigmaxx) &
      - suliciuvelocity*(righthsigmaxx-lefthsigmaxx)
    fluxright(5,:) = rightcells%hsigmazz*rightcells%velocity & 
      - rightvelocity*(rightcells%hsigmazz-righthsigmazz) &
      - suliciuvelocity*(righthsigmazz-lefthsigmazz)
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity>=0.) )
    fluxleft(1,:) = rightdepth*suliciuvelocity
    fluxleft(2,:) = rightdepth*suliciuvelocity**2 + suliciupressure
    fluxleft(3,:) = rightcells%tracer*rightdepth*suliciuvelocity
!     fluxleft(4,:) = leftcells%sigmaxx*leftcells%discharge &
!       + leftvelocity*(leftsigmaxx*leftdepth-leftcells%sigmaxx*leftcells%depth) &
!       + suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx)
!     fluxleft(5,:) = leftcells%sigmazz*leftcells%discharge &
!       + leftvelocity*(leftsigmazz*leftdepth-leftcells%sigmazz*leftcells%depth) &
!       + suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth)
    fluxleft(4,:) = leftcells%hsigmaxx*leftcells%velocity &
      + leftvelocity*(lefthsigmaxx-leftcells%hsigmaxx) &
      + suliciuvelocity*(righthsigmaxx-lefthsigmaxx)
    fluxleft(5,:) = leftcells%hsigmazz*leftcells%velocity &
      + leftvelocity*(lefthsigmazz-leftcells%hsigmazz) &
      + suliciuvelocity*(righthsigmazz-lefthsigmazz)
    fluxleft(6,:) = rightcells%microenstrophy*rightdepth*suliciuvelocity
    fluxleft(7,:) = rightcells%macroenstrophy*rightdepth*suliciuvelocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = rightcells%sigmaxx*rightcells%discharge & 
!       - rightvelocity*(rightcells%sigmaxx*rightcells%depth-rightsigmaxx*rightdepth)
!     fluxright(5,:) = rightcells%sigmazz*rightcells%discharge & 
!       - rightvelocity*(rightcells%sigmazz*rightcells%depth-rightsigmazz*rightdepth)
    fluxright(4,:) = rightcells%hsigmaxx*rightcells%velocity & 
      - rightvelocity*(rightcells%hsigmaxx-righthsigmaxx)
    fluxright(5,:) = rightcells%hsigmazz*rightcells%velocity & 
      - rightvelocity*(rightcells%hsigmazz-righthsigmazz)
    fluxright(6,:) = fluxleft(6,:)
    fluxright(7,:) = fluxleft(7,:)
  end where
  where( (suliciuvelocity<0.).AND.(rightvelocity<0.) )
    fluxleft(1,:) = rightcells%depth*rightcells%velocity
    fluxleft(2,:) = rightcells%depth*rightcells%velocity**2 + rightcells%pressure 
!     fluxleft(3,:) = rightcells%tracer*rightcells%discharge
!     fluxleft(4,:) = leftcells%sigmaxx*leftcells%discharge &
!       + leftvelocity*(leftsigmaxx*leftdepth-leftcells%sigmaxx*leftcells%depth) &
!       + suliciuvelocity*(rightsigmaxx*rightdepth-leftsigmaxx*leftdepth) &
!       + rightvelocity*(rightcells%sigmaxx*rightcells%depth-rightsigmaxx*rightdepth)
!     fluxleft(5,:) = leftcells%sigmazz*leftcells%discharge &
!       + leftvelocity*(leftsigmazz*leftdepth-leftcells%sigmazz*leftcells%depth) &
!       + suliciuvelocity*(rightsigmazz*rightdepth-leftsigmazz*leftdepth) &
!       + rightvelocity*(rightcells%sigmazz*rightcells%depth-rightsigmazz*rightdepth)
    fluxleft(3,:) = rightcells%htracer*rightcells%velocity
    fluxleft(4,:) = leftcells%hsigmaxx*leftcells%velocity &
      + leftvelocity*(lefthsigmaxx-leftcells%hsigmaxx) &
      + suliciuvelocity*(righthsigmaxx-lefthsigmaxx) &
      + rightvelocity*(rightcells%hsigmaxx-righthsigmaxx)
    fluxleft(5,:) = leftcells%hsigmazz*leftcells%velocity &
      + leftvelocity*(lefthsigmazz-leftcells%hsigmazz) &
      + suliciuvelocity*(righthsigmazz-lefthsigmazz) &
      + rightvelocity*(rightcells%hsigmazz-righthsigmazz)
    fluxleft(6,:) = rightcells%hmicroenstrophy*rightcells%velocity
    fluxleft(7,:) = rightcells%hmacroenstrophy*rightcells%velocity
    fluxright(1,:) = fluxleft(1,:)
    fluxright(2,:) = fluxleft(2,:)
    fluxright(3,:) = fluxleft(3,:)
!     fluxright(4,:) = rightcells%sigmaxx*rightcells%discharge
!     fluxright(5,:) = rightcells%sigmazz*rightcells%discharge
    fluxright(4,:) = rightcells%hsigmaxx*rightcells%velocity
    fluxright(5,:) = rightcells%hsigmazz*rightcells%velocity
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
! BELOW: parameterized derived data types would have been much useful (accessible only after fortran03)
  leftcells%depth = fluxleft(1,:)/leftcells%volume
  rightcells%depth = fluxright(1,:)/rightcells%volume
  leftcells%discharge = fluxleft(2,:)/leftcells%volume
  rightcells%discharge = fluxright(2,:)/rightcells%volume
  leftcells%htracer = fluxleft(3,:)/leftcells%volume
  rightcells%htracer = fluxright(3,:)/rightcells%volume
  leftcells%hsigmaxx = fluxleft(4,:)/leftcells%volume
  rightcells%hsigmaxx = fluxright(4,:)/rightcells%volume
  leftcells%hsigmazz = fluxleft(5,:)/leftcells%volume
  rightcells%hsigmazz = fluxright(5,:)/rightcells%volume
  leftcells%hmicroenstrophy = fluxleft(6,:)/leftcells%volume
  rightcells%hmicroenstrophy = fluxright(6,:)/rightcells%volume
  leftcells%hmacroenstrophy = fluxleft(7,:)/leftcells%volume
  rightcells%hmacroenstrophy = fluxright(7,:)/rightcells%volume
  where( rightvelocity<0 ) rightvelocity = -rightvelocity
  where( leftvelocity<0 ) leftvelocity = -leftvelocity
  t_neighbour = minval( (/rightcells%volume/rightvelocity,leftcells%volume/leftvelocity/) )
  ! ------------------------------------------------------------
  contains
  ! ------------------------------------------------------------
  subroutine suliciu_initialization( leftparameter, rightparameter, leftcells, rightcells, alphaspeed ) 
    use m_cell
    implicit none
    double precision, dimension(:), intent(inout) :: leftparameter, rightparameter  
    type(t_cell), dimension(:), intent(in) :: leftcells, rightcells
    double precision, intent(in) :: alphaspeed
    integer :: Nface
    double precision, dimension(:), allocatable :: temp0, temp1, temp2, temp3
    Nface = size(leftcells) !  = size(rightcells) = Ninterf
    allocate( temp0(Nface), temp1(Nface), temp2(Nface), temp3(Nface) )
    temp0 = rightparameter*rightcells%depth + leftparameter*leftcells%depth
    temp1 = rightcells%pressure-leftcells%pressure ! Pr-Pl
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
    temp1 = leftcells%velocity-rightcells%velocity ! u_l-u_r
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
