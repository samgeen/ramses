subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,emag_all,dt_all
  real(kind=8),dimension(4)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar+3),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  call velocity_fine(ilevel)

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  emag_all=0.0d0; emag_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        
        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do
        
        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,5)*vol
        end do
        
        ! Compute total magnetic energy
        do ivar=1,3
           do i=1,nleaf
              emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
        
        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,5)*vol
        end do
        do ivar=1,3
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
                   & -0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
        
        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if
        
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  comm_buffin(4)=emag_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,4,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc     ,dt_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
  emag_all=comm_buffout(4)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  emag_all=emag_loc
  dt_all=dt_loc
#endif

  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  emag_tot=emag_tot+emag_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
!  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::pi,time
  integer ::ivar,jgrid,ind_cell_vois
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,Cwnm
  real(dp)::dx_min, fact, Emag,Emag0


  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  Cwnm = sqrt(8000./scale_T2)

  pi=ACOS(-1.0d0)

  time = t * Cwnm / boxlen

 

  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  dx_min = (0.5D0**levelmin)/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        

       do i=1,ngrid


         if ( (xx(i,1) .lt. 4.*dx_min) .or. (boxlen-xx(i,1) .lt. 4.*dx_min) ) then 


           uold(ind_cell(i),1) = dens0

           uold(ind_cell(i),2) = 0.
           uold(ind_cell(i),3) = 0.
           uold(ind_cell(i),4) = 0.

           uold(ind_cell(i),6)  = sqrt(dens0*2.)*Cwnm*bx_bound
           uold(ind_cell(i),9)  = sqrt(dens0*2.)*Cwnm*bx_bound

           uold(ind_cell(i),7)  = sqrt(dens0*2.)*Cwnm*by_bound
           uold(ind_cell(i),10) = sqrt(dens0*2.)*Cwnm*by_bound

           uold(ind_cell(i),8)  = sqrt(dens0*2.)*Cwnm*bz_bound
           uold(ind_cell(i),11) = sqrt(dens0*2.)*Cwnm*bz_bound


          if( (xx(i,1) .gt. 2.*dx_min) .and. (xx(i,1) .lt. 4.*dx_min) ) then 

           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 

             !look for the grid neigbour of the left father
             jgrid = son(nbor(ind_grid(i),2))

             !remember iskip is calculated above
             !we must substract ngridmax because the neighbour of 4 (2) is 3 (1)
             ind_cell_vois = iskip + jgrid - ngridmax
           
             !this line is not useful since it is not used 
             ! it is there for consistency in the code
             uold(ind_cell(i),9) = uold(ind_cell_vois,6)


             !now calculate a factor to maintain the magnetic field reasonable

           A=0.5*(uold(ind_cell_vois,6)+uold(ind_cell_vois,nvar+1))
           B=0.5*(uold(ind_cell_vois,7)+uold(ind_cell_vois,nvar+2))
           C=0.5*(uold(ind_cell_vois,8)+uold(ind_cell_vois,nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)
           !Emag(t) / Emag(0)
           Emag0 = (dens0*8000./scale_T2*2.*1.5)*(bx_bound**2+by_bound**2+bz_bound**2)/2.
!           if (Emag0 .ne. 0) then
!               fact =  max(1., Emag / Emag0 - 1.)
!           else
               fact=1.
!           endif


             uold(ind_cell(i),2)=dens0*V0*Cwnm*(1.+ min(turb,1.)*(cos((xx(i,2)/boxlen-time)*5.*pi)*(cos((xx(i,3)/boxlen-time)*6.*pi+0.234)+0.2*cos((xx(i,2)/boxlen-time)*15.*pi))*cos((xx(i,3)/boxlen-time)*14.*pi+1.2)) ) / fact

             uold(ind_cell(i),3)=1.*dens0*turb*V0*Cwnm*( cos((xx(i,2)/boxlen-time)*7.*pi) * cos((xx(i,3)/boxlen-time)*8.*pi+2.15)  + 0.3*sin((xx(i,2)/boxlen-time)*17.*pi ) * sin((xx(i,3)/boxlen-time)*15.*pi+ 1.6 ) ) /fact 

             uold(ind_cell(i),4)=1.*dens0*turb*V0*Cwnm*( cos((xx(i,2)/boxlen-time)*6.*pi+1.23) * cos((xx(i,3)/boxlen-time)*5.*pi+1.76)  + 0.3*sin((xx(i,2)/boxlen-time)*16.*pi+2.1 ) * sin((xx(i,3)/boxlen-time)*17.*pi+2.645 ) )  /fact 

!((sin(xx(i,2)*pi/boxlen))**0.3)*((sin(xx(i,3)*pi/boxlen))**0.3)*

           endif

          endif


          if( (boxlen-xx(i,1) .gt. 2.*dx_min) .and. (boxlen-xx(i,1) .lt. 4.*dx_min) ) then 


           !because of the asymmetry introduced in umuscl (in ctoprim)  for the choice of the magnetic 
           ! component perpendicular to the face, we have to modify bx at the left face
           !ind should be either 1 or 3
          
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 

             !look for the grid neigbour of the left father
             jgrid = son(nbor(ind_grid(i),1))

             !remember iskip is calculated above
             !we must add ngridmax because the neighbour of 3 (1) is 4 (2)
             ind_cell_vois = iskip + jgrid + ngridmax
           
             uold(ind_cell(i),6) = uold(ind_cell_vois,9)


             !now calculate a factor to maintain the magnetic field reasonable

           A=0.5*(uold(ind_cell_vois,6)+uold(ind_cell_vois,nvar+1))
           B=0.5*(uold(ind_cell_vois,7)+uold(ind_cell_vois,nvar+2))
           C=0.5*(uold(ind_cell_vois,8)+uold(ind_cell_vois,nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

  
           ! Emag(t) / Emag(0)
           Emag0 = (dens0*8000./scale_T2*2.*1.5)*(bx_bound**2+by_bound**2+bz_bound**2)/2.
           if (Emag0 .ne. 0) then
               fact =  max(1., Emag / Emag0 - 1.)
           else
               fact=1.
           endif

             uold(ind_cell(i),2)=-dens0*V0*Cwnm*(1.+ min(turb,1.)*(cos((xx(i,2)/boxlen-time)*6.*pi+1.87)*(cos((xx(i,3)/boxlen-time)*5.*pi+0.78)+0.2*cos((xx(i,2)/boxlen-time)*14.*pi+0.23))*cos((xx(i,3)/boxlen-time)*13.*pi+1.67))) / fact

             uold(ind_cell(i),3)=1.*dens0*turb*V0*Cwnm*( cos((xx(i,2)/boxlen-time)*6.*pi+0.98) * cos((xx(i,3)/boxlen-time)*7.*pi+1.25)  + 0.3*sin((xx(i,2)/boxlen-time)*18.*pi+0.12 ) * sin((xx(i,3)/boxlen-time)*17.*pi+ 1.12 ) ) /fact 

             uold(ind_cell(i),4)=1.*dens0*turb*V0*Cwnm*( cos((xx(i,2)/boxlen-time)*5.*pi+0.23) * cos((xx(i,3)/boxlen-time)*6.*pi+0.7)  + 0.3*sin((xx(i,2)/boxlen-time)*17.*pi+1.01 ) * sin((xx(i,3)/boxlen-time)*16.*pi+2.1 ) ) /fact 


           endif

          endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           u=uold(ind_cell(i),2)/dens0
           v=uold(ind_cell(i),3)/dens0
           w=uold(ind_cell(i),4)/dens0

           !calculate total energy
           uold(ind_cell(i),neul)= dens0*Cwnm**2*gamma/(gamma-1.) + 0.5*dens0*(u**2+v**2+w**2) + Emag

         endif

       enddo



       
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
