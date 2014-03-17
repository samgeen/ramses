module cooling_module
  use cooling_module,only:X
  use rt_parameters

  implicit none
!=======================================================================
subroutine solve_cooling(nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none  
  integer::ncell
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost
  real(kind=8)::mu

  ! Neutral T = 100K
  ! Ionised T = 10000K
  mu = 1./(X*(1.+dU(2))

end subroutine solve_cooling
!=======================================================================

end module cooling_module
