!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp):: dmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,dd
  real(dp):: a1,a2,z0,a1_rho,a2_rho,G,mp,pi,Myear_s

  ! Some constants in cgs
  real(kind=8),parameter ::kB      = 1.3806200d-16
  !real(kind=8),parameter ::mH      = 1.6600000d-24
  !real(kind=8),parameter ::G       = 6.674d-8

!!  emass=gravity_params(1) ! Softening length
!!  xmass=gravity_params(2) ! Point mass coordinates
!!  ymass=gravity_params(3)
!!  zmass=gravity_params(4)
!!  dmass=1.0/(emass*(1.0+emass)**2)

  pi=ACOS(-1.)
  mp = 1.660531d-24  ! gramme
  G = 6.674d-8
  Myear_s=1.d6*365.*3600.*24.
  !Mpc=3.08567758d24
  Msolar = 1.9891d33

!! add the vertical galactic gravitational field
!! Kuijken & Gilmore 1989 taken from Joung & MacLow (2006)
  a1=1.42d-3 !!in kpc Myr-2 
  a2=5.49d-4 !!in Myr-2
  z0=0.18*1.d3 !!scale height in pc



!!the gravitational field is given by
!! g = -a1 z / sqrt(z^2+z0^2) - a2 z
!! rho = 1/(4piG) (a1 / z0) ( (z/z0)^2 + 1)^(-3/2) + a2/(4piG) 

  a1_rho = a1 / (4.*pi*G) / (z0/1.d3) / (Myear_s)**2 / mp

  a2_rho = a2 / (4.*pi*G)             / (Myear_s)**2 / mp

   ! Eqn 18 from Creasey et al 2013
   ! Calculate in cgs then convert
   ! Namelist parameters are fg, sigma_g
   T0=dx ! HACK - IMPORT T0 WITH dx THIS IS A BAD IDEA CHANGE THIS
   sigma_g_cgs = sigma_g * Msolar * pc ** (-2.0)
   invb_fact = mp*pi*G*sigma_g_cgs / (fg*kB*T0)

  do i=1,ncell
     rz=x(i,3)-0.5*boxlen
     d(i)= a1_rho / (1.+(rz/z0)**2)**(1.5) + a2_rho
     d(i) = 0.5 * sigma_g_cgs * invb_fact * &
            & cosh(rz*scale_l*invb_fact)**(-2.0) / scale_d
  end do

end subroutine rho_ana
