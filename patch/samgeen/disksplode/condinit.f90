!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_commons
  use amr_parameters
  use hydro_commons
  use poisson_parameters
!  use const
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  !TAKE CARE at this stage and  in this version this is not true for the
  ! first time that condinit is called
  ! because boxlen is determined in condinit
  ! for the first call x(i,1:3) are in  [0.,1.]
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer :: i,j,k,id,iu,iv,iw,ip
  real(dp):: pi
  integer :: ivar
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp),dimension(1:nvector)::rho     ! Calculate rho_ana with this

  real(dp),save:: first
  real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
  real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot
  integer,save:: n_size
  integer:: ind_i, ind_j, ind_k
  real(dp),save:: d_c,B_c,ind,seed1,seed2,seed3,xi,yi,zi,zeta
  real(dp),save:: res_int,r_0,C_s,omega,v_rms,cont_ic,mass_total,mass_tot2,min_col_d,max_col_d
  real(dp):: col_d,eli,sph,vx,vy,vz
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::xl,yl,zl,xx,yy,zz,bx,by,bz,dxmin
  integer:: ii,jj,kk,nticks
  real(dp)::ener_rot,ener_grav,ener_therm,ener_grav2,ener_turb
  real(dp),dimension(1000):: mass_rad
  real(dp):: mu=1.4d0 ! NOTE - MUST BE THE SAME AS IN units.f90!!
  real(dp):: T0=1e4
!  real(dp)::myid
  real(dp)::P_WNM=0.0d0


!    myid=1

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu


   !do various things which needs to be done only one time
   if( first .eq. 0.) then
    id=1; iu=2; iv=3; iw=4; ip=5
    pi=acos(-1.0d0)



    if(myid==1) write(*,*) '** ENTER  in condinit **'


    call calc_boxlen

#ifdef TURBULENCE
    !now read the turbulent velocity field used as initial condition
    if( myid ==1) write(*,*) 'Read the file which contains the initial turbulent velocity field'
    open(20,file='ramses.data',form='formatted')
    read(20,*) n_size, ind, seed1,seed2,seed3

     if(n_size .ne. 100) then
       write(*,*) 'Unextected field size'
       stop
     endif

     v_rms=0.
     mass_total=0.
     mass_tot2 =0.
     do k=1,n_size
     do j=1,n_size
     do i=1,n_size
        read(20,*)xi,yi,zi,vx,vy,vz
        q_idl(1,i,j,k) = vx
        q_idl(2,i,j,k) = vy
        q_idl(3,i,j,k) = vz

        xi = boxlen/bl_fac*((i-0.5)/n_size-0.5)
        yi = boxlen/bl_fac*((j-0.5)/n_size-0.5)
        zi = boxlen/bl_fac*((k-0.5)/n_size-0.5)
        eli =  (xi/r_0)**2+(yi/r_0)**2+(zi/(r_0*rap))**2
        if( eli .lt. zeta**2) then

         vx_tot = vx_tot + d_c/(1.+eli)*vx
         vy_tot = vy_tot + d_c/(1.+eli)*vy
         vz_tot = vz_tot + d_c/(1.+eli)*vz

         vx2_tot = vx2_tot + d_c/(1.+eli)*vx**2
         vy2_tot = vy2_tot + d_c/(1.+eli)*vy**2
         vz2_tot = vz2_tot + d_c/(1.+eli)*vz**2

         ener_turb = ener_turb + d_c/(1.+eli)*(vx**2+vy**2+vz**2)
         mass_total = mass_total +  d_c / (1.+eli)
         ener_rot = ener_rot + d_c/(1.+eli) * omega**2 * (yi**2 + zi**2)
        endif
     enddo
!       eli = (yi/r_0)**2 + (zi/r_0/rap)**2
!        if( eli .lt. zeta**2) then
!          col_d = r_0*d_c/sqrt(1.+eli)*atan( sqrt( (zeta**2-eli)/(1.+eli) ) )
!          mass_tot2 = mass_tot2 + col_d
!        endif
     enddo
     enddo
    close(20)

     vx_tot = vx_tot / mass_total
     vy_tot = vy_tot / mass_total
     vz_tot = vz_tot / mass_total

     vx2_tot = vx2_tot / mass_total
     vy2_tot = vy2_tot / mass_total
     vz2_tot = vz2_tot / mass_total

     v_rms = sqrt( vx2_tot-vx_tot**2 + vy2_tot-vy_tot**2 + vz2_tot-vz_tot**2 )

     mass_total = mass_total*(boxlen/n_size)**3
     if (myid ==1) write(*,*) 'We verify the calculation for the mass. The 2 following values must be very close:'
     if (myid ==1) write(*,*) 'mass_total, mass_c ',mass_total, mass_c !,mass_tot2

     ener_rot  = 0.5 * ener_rot*(boxlen/n_size)**3
     ener_turb = 0.5 * ener_turb*(boxlen/n_size)**3

     !estimate of the thermal over gravitational energy
     if (myid == 1) write(*,*) 'estimate (uniform density is assumed) of the ratio of thermal over gravitational energy'
     if (myid == 1) write(*,*)  ener_therm / ener_grav

     if (myid == 1) write(*,*) 'good estimate of the ratio of thermal over gravitational energy'
     if (myid == 1) write(*,*)  ener_therm / ener_grav2

     !estimate of the rotational over gravitational energy ratio
     if (myid .eq. 1) write(*,*) 'estimate of the rotational over gravitational energy ratio'
     if (myid .eq. 1) write(*,*) 'ener_rot/ener_grav2 ', ener_rot / ener_grav2


     !calculate now the coefficient by which the turbulence velocity needs
     !to be multiplied

     if (myid .eq. 1) write(*,*) 'vrms non norm ',v_rms

    !ph 01/09 new definition entails r_0 instead of r_0 * zeta, the external radius
     v_rms = ff_vct * sqrt(32.*d_c/3./pi)*r_0 / v_rms

     if (myid .eq. 1) write(*,*) 'vrms mult ',v_rms


     !estimate of the turbulent over gravitational energy ratio
     if (myid .eq. 1) write(*,*) 'estimate of the turbulent over gravitational energy ratio'
     if (myid .eq. 1) write(*,*) 'ener_turb/ener_grav2 ', ener_turb*(v_rms**2) / ener_grav2



    100 format(i5,4e10.5)
    101 format(6e10.5)
    102 format(i5)

    if (myid ==1)  write(*,*) 'Reading achieved'

#endif
    first = 1.
   endif


   ! Get the density field
   ! HACK - pass T0 to this by replacing dx
   T0 = 1e4 ! TODO: change cooling function
   call rho_ana(x,rho,T0,nn)
   


   DO i=1,nn


       x(i,1) = x(i,1) - 0.5*boxlen
       x(i,2) = x(i,2) - 0.5*boxlen
       x(i,3) = x(i,3) - 0.5*boxlen


       !initialise the density field
       ! Eqn 16 in Creasey et al 2013
       q(i,1) = rho(i)
       q(i,2:4) = 0.0
       q(i,5) = q(i,1) * T0 / scale_T2
       q(i,6:8) = 0.0

#ifdef TURBULENT

       if(all(abs(x(i,:) / (boxlen / bl_fac)) <= 0.5)) then
         !initialise the turbulent velocity field
         !make a zero order interpolation (should be improved)
         ind_i = int((x(i,1)/(boxlen/bl_fac)+0.5)*n_size)+1
         ind_j = int((x(i,2)/(boxlen/bl_fac)+0.5)*n_size)+1
         ind_k = int((x(i,3)/(boxlen/bl_fac)+0.5)*n_size)+1


         if( ind_i .lt. 1 .or. ind_i .gt. n_size) write(*,*) 'ind_i ',ind_i,boxlen,x(i,1),n_size
         if( ind_j .lt. 1 .or. ind_j .gt. n_size) write(*,*) 'ind_j ',ind_j
         if( ind_k .lt. 1 .or. ind_k .gt. n_size) write(*,*) 'ind_k ',ind_k

         ! STG HACK - remove vel
         !q(i,2) = v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot)
         !q(i,3) = v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot)
         !q(i,4) = v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot)

         !add  rotation. x cos(thet_mag) + y sin(thet_mag) is the rotation axis
         sph =  (x(i,1)/r_0)**2+(x(i,2)/r_0)**2+(x(i,3)/(r_0))**2
         if( sph .lt. (zeta*rap)**2 ) then

         !to check these formulae one can verify that those arrays are perpendicular
         !with (cos(thet_mag),sin(thet_mag),0) and that the norm of the vectorial product of
         ! (cos(thet_mag),sin(thet_mag),0) by the above arrays is equal to the distance
         !  (x sin(thet)-y cos(thet))^2 + z^2
         ! STG HACK - remove vel
         !  q(i,2) = q(i,2) - (omega*x(i,3)*sin(thet_mag))
         !  q(i,3) = q(i,3) + (omega*x(i,3)*cos(thet_mag))
         !  q(i,4) = q(i,4) + (omega*(x(i,1)*sin(thet_mag)-x(i,2)*cos(thet_mag)))
         endif
       else
          ! STG HACK - remove vel
         !q(i, 2:4) = 0.0
       endif

#endif


  ENDDO


  dxmin=boxlen*0.5d0**(levelmin+3)

  if( dx .lt. dxmin) then
    write(*,*) 'dxmin too large'
    write(*,*) 'dx ',dx/boxlen
    write(*,*) 'dxmin ',dxmin/boxlen
    stop
  endif

  nticks=dx/dxmin


! STG HACK - remove MHD field
#ifdef SOLVERmhdHACK
   DO i=1,nn
   q(i,6)=0.

     xl=x(i,1)-0.5*dx
     yl=x(i,2)-0.5*dx
     zl=x(i,3)-0.5*dx

     !the magnetic field in cells must be subdivided in order to insure that the magnetic
     !flux is the same in coarse and refined grids
     DO jj=1,nticks
     DO kk=1,nticks

        yy=yl+(dble(jj)-0.5d0)*dxmin
        zz=zl+(dble(kk)-0.5d0)*dxmin

       !this formula comes from the integration of the density distribution along x
       eli = (yy/r_0)**2 + (zz/r_0/rap)**2
        if( eli .lt. zeta**2) then
         col_d = r_0*d_c/sqrt(1.+eli)*atan( sqrt( (zeta**2-eli)/(1.+eli) ) )
         col_d = max(col_d,min_col_d)
        else
         col_d = min_col_d
        endif

       !Bx component
       q(i,6     ) = q(i,6) + B_c * col_d / max_col_d
       q(i,nvar+1) = q(i,6)

       !By component
       q(i,7     ) = 0.
       q(i,nvar+2) = 0.

       !Bz component
       q(i,8     ) = 0.
       q(i,nvar+3) = 0.

     ENDDO
     ENDDO

       q(i,6:8)           = q(i,6:8)           / dble(nticks)**2

!new version rotates the rotation velocity
     !rotates the magnetic field of an angle theta
!       bx=q(i,6)
!       by=q(i,7)
!       q(i,6) =  bx*cos(thet_mag) + by*sin(thet_mag)
!       q(i,7) =  bx*sin(thet_mag) - by*cos(thet_mag)

       q(i,nvar+1:nvar+3) = q(i,6:8)
  ENDDO
#endif

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  !u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  !u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  !u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  !kinetic + magnetic energy
!#ifdef SOLVERmhd
!  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,6)+q(1:nn,nvar+1))**2
!  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,7)+q(1:nn,nvar+2))**2
!  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,8)+q(1:nn,nvar+3))**2
!#endif
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic field
#ifdef SOLVERmhd
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#else
  ! passive scalars
  do ivar=6,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================

subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)

     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana
!================================================================
!================================================================
!================================================================
!================================================================
subroutine calc_dmin(d_c)
  use amr_commons
  use hydro_commons
  implicit none

  real(dp):: d_c, cont_ic

  cont_ic = 10.
  dmin = d_c / cont / cont_ic

  if (myid == 1) then
    write(*,*) "dmin = ", dmin
  endif
end subroutine calc_dmin
!================================================================
!================================================================
!================================================================
!================================================================
subroutine calc_boxlen
  use amr_commons
  use amr_parameters
  use hydro_commons
  use poisson_parameters
!  use const
  implicit none
  !================================================================
  !this routine calculate boxlen
  !================================================================
  integer :: i
  real(dp):: pi
  real(dp):: d_c,zeta
  real(dp):: res_int,r_0,C_s,kpc
  integer::  np
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),save:: first
  real(dp):: mu=1.4d0 ! NOTE - MUST BE THE SAME AS IN units.f90!!
!  real(dp)::myid

!   myid=1

    if (first .eq. 0.) then

    pi=acos(-1.0d0)

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    ! Set boxlen to 1kpc
    kpc=3.08567758d21
    boxlen = kpc/scale_l


    first=1.
    endif

    call calc_dmin(d_c)

end subroutine calc_boxlen
