module m_cell
  implicit none
  type t_cell
    double precision center, volume
    double precision depth
    double precision discharge, velocity
    double precision tracer, sigmaxx, sigmazz
    double precision pressure
  end type t_cell
  ! creer un operateur << ?
end module m_cell
