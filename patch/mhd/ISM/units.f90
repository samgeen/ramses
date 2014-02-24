subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module

  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l,mp,G,pc,mu,kbol
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------
  mu  =  1.4
  mp =  mu * 1.660531d-24  ! gramme
  G = 6.67d-8
  kbol  =  1.38062d-16   ! erg/degre

  ! scale_d converts mass density from user units into g/cc
  scale_d = mp  !code density units is in particle/cc

  ! scale_t converts time from user units into seconds
  scale_t = 1.0/sqrt(G*scale_d)

  ! scale_l converts distance from user units into cm
  scale_l = 3.08d18 !code distance units is in parsec

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into T in Kelvin
  scale_T2 = mp/kbol * scale_v**2

!   write(*,*) scale_d,scale_t,scale_l,scale_v,scale_T2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mp * scale_d

end subroutine units
