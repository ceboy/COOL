module m_cell
  implicit none
  type t_cell
    double precision center, volume
    double precision depth, discharge, velocity, pressure, tracer
    double precision sigmaxx, sigmazz
  end type t_cell
  ! creer un operateur << ?
end module m_cell
