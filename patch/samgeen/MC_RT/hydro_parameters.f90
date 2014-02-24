module hydro_parameters
  use amr_parameters

  ! Number of independant variables
#ifndef NVAR
  integer,parameter::nvar=ndim+2
#else
  integer,parameter::nvar=NVAR
#endif
  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
#if NVAR > NDIM+2
  real(dp),dimension(1:MAXBOUND,1:NVAR-NDIM-2)::var_bound=0.0
#endif
  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
#if NVAR > NDIM+2
  real(dp),dimension(1:NVAR-NDIM-2)::err_grad_var=-1.0
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
#if NVAR > NDIM+2
  real(dp),dimension(1:MAXREGION,1:NVAR-NDIM-2)::var_region=0.0
#endif

  !initial temperature used for the isothermal run
  real(dp)::temper
  real(dp)::temper_iso

  !Initial conditions parameter for the dense core
  real(dp)::bx_bound=0.
  real(dp)::by_bound=0.
  real(dp)::bz_bound=0.
  real(dp)::turb=0.
  real(dp)::dens0=0.
  real(dp)::V0=0.
  real(dp)::Height0=0.

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::gamma=1.4d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Passive variables index
  integer::imetal=6
  integer::idelay=6
  integer::ixion=6
  integer::ichem=6

  ! Supernova list
  integer, parameter:: MAXSN=100
  integer:: sn_count, sn_npart, sn_i = 1
  logical:: sn_all = .false.
  real(dp), dimension(1:MAXSN):: sn_radius=1.0, sn_mass=0.0, sn_energy=0.0, sn_time=0.0, sn_part_radius=1.0
  real(dp), dimension(1:MAXSN, 3):: sn_center=0.0
  real(dp):: sn_freq=0.
  real(dp):: sn_freq_mult=0.
  real(dp)::t_last_sn
  real(dp)::sn_e_ref
  real(dp)::sn_mass_ref
  logical::supernovae=.false.

  !Initial conditions parameter for the dense core
  real(dp)::mass_c=1.   !cloud mass in solar mass
  real(dp)::rap=1.      !axis ratio
  real(dp)::cont=1.     !density contras
  real(dp)::ff_sct=1.   !freefall time / sound crossing time
  real(dp)::ff_rt=1.    !freefall time / rotation time
  real(dp)::ff_act=1.   !freefall time / Alfven crossing time
  real(dp)::ff_vct=1.   !freefall time / Vrms crossing time
  real(dp)::thet_mag=0. !angle between magnetic field and rotation axis
  real(dp)::bl_fac=1.   !multiply calculated boxlen by this factor


end module hydro_parameters
