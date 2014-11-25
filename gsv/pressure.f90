subroutine pressure(cell,theta,gavrilyuk,elasticmodulus,oneoverell)
  use m_physics
  use m_cell
  implicit none
  type (t_cell), dimension(:), intent(inout) :: cell 
  double precision, intent(in) :: theta, gavrilyuk, elasticmodulus, oneoverell
!   integer :: Nx
!   Nx = size(cell)-2
!   cell(2:Nx+1)%pressure = cos(theta)*g*cell(2:Nx+1)%depth**2/2 &
!     + gavrilyuk*cell(2:Nx+1)%depth**3*( cell(2:Nx+1)%microenstrophy + cell(2:Nx+1)%macroenstrophy ) &
!     + elasticmodulus*(cell(2:Nx+1)%sigmazz-cell(2:Nx+1)%sigmaxx)*cell(2:Nx+1)%depth/ &
! 	(1. + oneoverell*(cell(2:Nx+1)%sigmazz+cell(2:Nx+1)%sigmaxx))
!   cell(2:Nx+1)%speed = sqrt( cos(theta)*g*cell(2:Nx+1)%depth &
!     + gavrilyuk*3*cell(2:Nx+1)%depth**2*( cell(2:Nx+1)%microenstrophy + cell(2:Nx+1)%macroenstrophy ) &
!     + elasticmodulus*(3*cell(2:Nx+1)%sigmazz+cell(2:Nx+1)%sigmaxx)/ &
! 	(1. - oneoverell*(cell(2:Nx+1)%sigmazz+cell(2:Nx+1)%sigmaxx)) &
!     + oneoverell*2*elasticmodulus*((cell(2:Nx+1)%sigmazz-cell(2:Nx+1)%sigmaxx)/ &
! 	(1. - oneoverell*(cell(2:Nx+1)%sigmazz+cell(2:Nx+1)%sigmaxx)))**2 )
  cell%pressure = cos(theta)*g*cell%depth**2/2 &
    + gavrilyuk*cell%depth**3*( cell%microenstrophy + cell%macroenstrophy ) &
    + elasticmodulus*(cell%sigmazz-cell%sigmaxx)*cell%depth/ &
	(1. + oneoverell*(cell%sigmazz+cell%sigmaxx))
  cell%speed = sqrt( cos(theta)*g*cell%depth &
    + gavrilyuk*3*cell%depth**2*( cell%microenstrophy + cell%macroenstrophy ) &
    + elasticmodulus*(3*cell%sigmazz+cell%sigmaxx)/ &
	(1. - oneoverell*(cell%sigmazz+cell%sigmaxx)) &
    + oneoverell*2*elasticmodulus*((cell%sigmazz-cell%sigmaxx)/ &
	(1. - oneoverell*(cell%sigmazz+cell%sigmaxx)))**2 )
end subroutine pressure
