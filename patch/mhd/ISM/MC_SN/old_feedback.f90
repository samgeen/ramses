subroutine make_sn
  use amr_commons
  use hydro_commons
  implicit none

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx
  real(dp):: scale, dx_min, dx_loc, vol_loc
  real(dp), dimension(1:3):: xbound, skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok

  real(dp),dimension(1:ndim):: sn_cent
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_e, sn_vol, sn_d, sn_ed, sn_rp, dx_sel
  real(dp):: rr, pi
  logical:: sel = .false.
  logical, save:: first = .true.
  real(dp)::xseed

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_sn'

  pi = acos(-1.0)

  ! Mesh spacing in that level
  xbound(1:3) = (/ dble(nx), dble(ny), dble(nz) /)
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
  skip_loc(1) = dble(icoarse_min)
  skip_loc(2) = dble(jcoarse_min)
  skip_loc(3) = dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  if (first) then 
     xseed=0.5
     call random_number(xseed)
     first=.false.
  endif


  if(sn_freq_mult .eq. 0.) then
     sn_r = sn_radius(sn_i)
     sn_m = sn_mass(sn_i)
     sn_e = sn_energy(sn_i)
     sn_rp = sn_part_radius(sn_i)
     sn_cent(1) = sn_center(sn_i, 1)
     sn_cent(2) = sn_center(sn_i, 2)
     sn_cent(3) = sn_center(sn_i, 3)
  else
     sn_r = 3.*(0.5**levelmin)*scale
     sn_m = sn_mass_ref
     sn_e = sn_e_ref 
     sn_rp = 0.
     call random_number(xseed)
     sn_cent(1)= xseed*boxlen
     call random_number(xseed)
     sn_cent(2)= xseed*boxlen
     call random_number(xseed)
     sn_cent(3)=(0.5 + (xseed-0.5) / 10. )*boxlen !assume mid-plane for the moment 
  endif

  if(sn_r /= 0.0) then
    sn_vol = 4. / 3. * pi * sn_r**3
    sn_d = sn_m / sn_vol
    sn_ed = sn_e / sn_vol
  end if

  if(myid .eq. 1 .or. myid .eq. 2) then
   write(*,*) 'x_sn, y_sn, z_sn, ',sn_cent(1),sn_cent(2),sn_cent(3)
  endif

  ! Loop over levels
  do ilevel = levelmin, nlevelmax
    ! Computing local volume (important for averaging hydro quantities)
    dx = 0.5d0**ilevel
    dx_loc = dx * scale
    vol_loc = dx_loc**ndim
    !if(.not. sel) then
      ! dx_sel will contain the size of the biggest leaf cell around the center
      !dx_sel = dx_loc
      !sn_vol = vol_loc
    !end if

    ! Cell center position relative to grid center position
    do ind=1,twotondim
      iz = (ind - 1) / 4
      iy = (ind - 1 - 4 * iz) / 2
      ix = (ind - 1 - 2 * iy - 4 * iz)
      xc(ind,1) = (dble(ix) - 0.5d0) * dx
      xc(ind,2) = (dble(iy) - 0.5d0) * dx
      xc(ind,3) = (dble(iz) - 0.5d0) * dx
    end do

    ! Loop over grids
    ncache=active(ilevel)%ngrid
    do igrid = 1, ncache, nvector
      ngrid = min(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

      ! Loop over cells
      do ind = 1, twotondim
        ! Gather cell indices
        iskip = ncoarse + (ind - 1) * ngridmax
        do i = 1, ngrid
          ind_cell(i) = iskip + ind_grid(i)
        end do

        ! Gather cell center positions
        do i = 1, ngrid
          xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
        end do
        ! Rescale position from coarse grid units to code units
        do i = 1, ngrid
           xx(i, :) = (xx(i, :) - skip_loc(:)) * scale
        end do

        ! Flag leaf cells
        do i = 1, ngrid
          ok(i) = (son(ind_cell(i)) == 0)
        end do

        do i = 1, ngrid
          if(ok(i)) then
            if(sn_r == 0.0) then
              sn_d = sn_m / vol_loc ! XXX
              sn_ed = sn_e / vol_loc ! XXX
              rr = 1.0
              do idim = 1, 3
                !rr = rr * max(1.0 - abs(xx(i, idim) - sn_center(sn_i, idim)) / dx_sel, 0.0)
                rr = rr * max(1.0 - abs(xx(i, idim) - sn_cent(idim)) / dx_loc, 0.0)
              end do
              !if(rr > 0.0) then
                !if(.not. sel) then
                  !! We found a leaf cell near the supernova center
                  !sel = .true.
                  !sn_d = sn_m / sn_vol
                  !sn_ed = sn_e / sn_vol
                !end if
                uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d * rr
                uold(ind_cell(i), 5) = uold(ind_cell(i), 5) + sn_ed * rr
              !end if
            else
              rr = sum(((xx(i, :) - sn_cent(:)) / sn_r)**2)

              if(rr < 1.) then
                uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                uold(ind_cell(i), 5) = uold(ind_cell(i), 5) + sn_ed
              endif
            endif
          endif
        end do
      end do
      ! End loop over cells
    end do
    ! End loop over grids
  end do
  ! End loop over levels

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo
end subroutine make_sn
