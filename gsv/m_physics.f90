module m_physics
  implicit none
  !--------------------------------------------------------------------------
  double precision, parameter :: Pi = 3.1415927
  double precision, parameter :: g = 9.81d0 ! gravity constant
  double precision :: k = 0d0  ! friction parameter
  double precision :: kwave = 2.5  ! wave number
end module m_physics