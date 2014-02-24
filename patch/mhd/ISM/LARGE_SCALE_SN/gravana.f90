!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz
  real(dp)::a1,a2,z0,z02,Myear_s,z,omega,scale
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,nx_loc

  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then 
     gmass=gravity_params(1) ! GM
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
        f(i,1)=-gmass*rx/rr**3
#if NDIM>1
        f(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        f(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if


  ! galactic field
  if(gravity_type==3)then 

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
!  scale_kappa=1./scale_l
!  scale_m=scale_d*scale_l**3

  !compute galactic gravitational acceleration
  Myear_s=1.d6*365.*3600.*24.

  nx_loc=(icoarse_max-icoarse_min+1)

  scale=boxlen/dble(nx_loc)

!! add the vertical galactic gravitational field
!! Kuijken & Gilmore 1989 taken from Joung & MacLow (2006)
  a1=1.42d-3 !!in kpc Myr-2 
  a2=5.49d-4 !!in Myr-2
  z0=0.18*1.d3 !!scale height in pc

!!now convert in code units assuming boxlen is in pc
  a1=a1*1.d3/(Myear_s**2)*(scale_t**2)
  a2=a2/(Myear_s**2)*(scale_t**2)
  z02=z0**2  




     do i=1,ncell

       !!!add here the galactic acceleration on the gas
       !!the gravitational field is given by
       !! g = -a1 z / sqrt(z^2+z0^2) - a2 z
       z = x(i,3)-0.5*scale
       
       f(i,1) = 0.
       f(i,2) = 0.

       f(i,3) =  - (a1/sqrt(z**2+z02)+a2)*z
     end do
  end if


end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
