!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
!  use hydro_parameters
  use hydro_commons
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i,j,k
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::pi,xx,yy
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,mag_norm

  real(dp),save:: first
  real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
  real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot,vx,vy,vz,v_rms
  integer,save:: n_size
  integer:: ind_i, ind_j, ind_k
  real(dp),save:: ind,seed1,seed2,seed3,xi,yi,zi
  real(dp):: n_total


  ! Call built-in initial condition generator
!  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........
  pi=ACOS(-1.0d0)


  mass_sph = 10. * (boxlen*(0.5**levelmin))**3

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
!  write(*,*) 'Echelle de temperature',scale_T2
!  write(*,*) 'temperature du code',8000./scale_T2
!  write(*,*) 'dens0 ',dens0

   mag_norm = sqrt(1.*8000./scale_T2*2.*1.5)

   if( first .eq. 0.) then 

      temper = (8000. / scale_T2 ) / dens0

     if (isothermal) temper = temper_iso / scale_T2


    !now read the turbulent velocity field used as initial condition
    if( myid ==1) write(*,*) 'Read the file which contains the initial turbulent velocity field'
    open(20,file='ramses.data',form='formatted')
    read(20,*) n_size, ind, seed1,seed2,seed3

     if(n_size .ne. 100) then 
       write(*,*) 'Unextected field size'
       stop
     endif
     
     v_rms=0.

     vx_tot  = 0.
     vy_tot  = 0.
     vz_tot  = 0.
     vx2_tot = 0.
     vy2_tot = 0.
     vz2_tot = 0.

     do k=1,n_size
     do j=1,n_size
     do i=1,n_size
        read(20,*)xi,yi,zi,vx,vy,vz
        q_idl(1,i,j,k) = vx
        q_idl(2,i,j,k) = vy
        q_idl(3,i,j,k) = vz

        xi = boxlen*((i-0.5)/n_size-0.5)
        yi = boxlen*((j-0.5)/n_size-0.5)
        zi = boxlen*((k-0.5)/n_size-0.5)

         vx_tot = vx_tot + vx
         vy_tot = vy_tot + vy
         vz_tot = vz_tot + vz

         vx2_tot = vx2_tot + vx**2
         vy2_tot = vy2_tot + vy**2
         vz2_tot = vz2_tot + vz**2

     enddo
     enddo
     enddo
    close(20)

    n_total = n_size**3

     vx_tot = vx_tot / n_total
     vy_tot = vy_tot / n_total
     vz_tot = vz_tot / n_total

     vx2_tot = vx2_tot / n_total
     vy2_tot = vy2_tot / n_total
     vz2_tot = vz2_tot / n_total

     v_rms = sqrt( vx2_tot-vx_tot**2 + vy2_tot-vy_tot**2 + vz2_tot-vz_tot**2 ) 

     

     !calculate now the coefficient by which the turbulence velocity needs
     !to be multiplied 


     !turb is in km/s ,  1.d5 converts it in cm/s
     v_rms =  turb*1.d5 / scale_v / v_rms

    if (myid ==1 ) write(*,*), 'turb ', turb, ', v_rms ', v_rms , 'first ',first


    100 format(i5,4e10.5)
    101 format(6e10.5)
    102 format(i5)

    if (myid ==1)  write(*,*) 'Reading achieved'
    first = 1.
   endif


!     write(*,*), 'turb ', turb, ', v_rms ', v_rms , 'dens_0, mag_norm ',dens0,mag_norm

  do i=1,nn

     x(i,1) = x(i,1) - 0.5*boxlen
     x(i,2) = x(i,2) - 0.5*boxlen
     x(i,3) = x(i,3) - 0.5*boxlen

     !Bx component
     q(i,6     ) = bx_bound * mag_norm * exp(-x(i,3)**2/(2.*Height0**2)) !exponential profile along z

     q(i,nvar+1) = q(i,6)

     !By component
     q(i,7     ) = 0. !by_bound * mag_norm
     q(i,nvar+2) = 0. !by_bound * mag_norm

     !Bz component
     q(i,8     ) = 0. !bz_bound * mag_norm
     q(i,nvar+3) = 0. !bz_bound * mag_norm

     !en cgs
        !densite
     q(i,1) = dens0 * max(exp(-x(i,3)**2/(2.*Height0**2)),1.d-2) !exponential profile along z
        !pression
     q(i,5) =  q(i,1)*temper !(8000. / scale_T2 )

!     temper = q(i,5) / q(i,1)


       !initialise the turbulent velocity field
       !make a zero order interpolation (should be improved)
       ind_i = int((x(i,1)/boxlen+0.5)*n_size)+1
       ind_j = int((x(i,2)/boxlen+0.5)*n_size)+1
       ind_k = int((x(i,3)/boxlen+0.5)*n_size)+1


       if( ind_i .lt. 1 .or. ind_i .gt. n_size) write(*,*) 'ind_i ',ind_i,boxlen,x(i,1),n_size
       if( ind_j .lt. 1 .or. ind_j .gt. n_size) write(*,*) 'ind_j ',ind_j
       if( ind_k .lt. 1 .or. ind_k .gt. n_size) write(*,*) 'ind_k ',ind_k

       q(i,2) = v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot)
       q(i,3) = v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot)
       q(i,4) = v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot)


  end do


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  !kinetic + magnetic energy
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,8)+q(1:nn,nvar+3))**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic field 
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do


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
