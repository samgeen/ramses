subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(active(ilevel)%ngrid==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if(cooling.and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,neul=5
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,Temp,dt_ilev,divB
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekin,emag,T2min,Zsolar

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=abs(uold(ind_leaf(i),1))
     end do
     
     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),neul)
     end do
     do i=1,nleaf
        ekin(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           ekin(i)=ekin(i)+0.5d0*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        emag(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+neul)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do

     ! Compute temperature 
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekin(i)-emag(i))/nH(i)
     end do

     if (isothermal) then 

       !impose isothermallity
     do i=1,nleaf
!       write(*,*) 'T2min',T2min(i),temper_iso, scale_T2
       T2min(i) = temper_iso / scale_T2
        if( nH(i) .le. smallr) then
          ekin(i) = 0.
          do idim=1,3
            uold(ind_leaf(i),idim+1) = uold(ind_leaf(i),idim+1) * smallr / abs(uold(ind_leaf(i),1))
            uold(ind_leaf(i),1) = smallr
            ekin(i)=ekin(i)+0.5d0*uold(ind_leaf(i),idim+1)**2 / smallr 
          enddo
        endif
     end do

     else

     ! Update temperature from cooling function
     do i=1,nleaf
        Temp = T2(i)
        dt_ilev = dtnew(ilevel)
        if( Temp .lt. 0.) then 
!         write(*,*) 'NN, TT ',nH(i), Temp
!         write(*,*) uold(ind_leaf(i),2),uold(ind_leaf(i),3)
!         write(*,*) 'ekin, emag, etot ',ekin(i),emag(i),uold(ind_leaf(i),5)
        endif
        call calc_temp(nH(i),Temp,dt_ilev)
        T2min(i) = Temp
        if( nH(i) .eq. smallr) then
          ekin(i) = 0.
          do idim=1,3
            uold(ind_leaf(i),idim+1) = uold(ind_leaf(i),idim+1) * smallr / abs(uold(ind_leaf(i),1))
            uold(ind_leaf(i),1) = nH(i)
            ekin(i)=ekin(i)+0.5d0*uold(ind_leaf(i),idim+1)**2/nH(i)
          enddo
        endif
!       if ( mod(i,100) .eq. 0) write(*,*) i,T2(i),T2min(i)
     end do

     endif


     ! compute total energy
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/(gamma-1.0) + ekin(i) + emag(i)
     end do

     ! Update total fluid energy
     do i=1,nleaf
        uold(ind_leaf(i),neul) = T2min(i)
     end do

  end do
  ! End loop over cells

end subroutine coolfine1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  calc_temp(NN,TT,dt_tot_unicode)
  use amr_parameters
  use hydro_commons

  implicit none

  integer :: n,i,j,k,idim, iter, itermax,ii

  real(dp) :: dt, dt_tot, temps, dt_max, itermoy
  real(dp) :: rho,temp,dt_tot_unicode

  real(dp) :: mm,uma, kb, alpha,mu,kb_mm
  real(dp) :: NN,TT, TTold, ref,dRefdT, eps, vardt,varrel, dTemp
  real(dp) :: rhoutot2
  real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
!
! Cette routine fonctionne en cgs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   kb  =  1.38062d-16   ! erg/degre
!  uma =  1.660531e-24  ! gramme
!  mu  =  1.4
!  mm = mu*uma 
!  kb_mm = kb / mm 
!  TT = TT  / kb  !/ kb_mm 


  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)


   if( TT .le. 0.) then 
     TT = 50. / scale_T2
     write(*,*) 'negative temperature ',TT
     return
   endif  


!   !use for the shell problem to keep high the temperature of the coronal gas 
!   if( TT * scale_T2 .ge. 20000. ) then 
!     return
!   endif  



!if( TT*scale_T2 .gt. 50.) then
!TT = 50. / scale_T2
!return
!endif


  vardt = 10.**(1./10.); varrel = 0.2

 

   dt_tot = dt_tot_unicode * scale_t ! * 3.08d18 / sqrt(kb_mm)
   TT     = TT * scale_T2


!  nn = (rho/(gramme/cm3)) /mm
  
  itermax = 0 ; itermoy = 0.



        if (NN .le. smallr) then 
         if( NN .le. 0)  write(*,*) 'prob dens',NN
         NN = smallr  !max(NN,smallr)
        endif


!     alpha = NN*kb_mm/(gamma-1.)
     alpha = NN*kb/(gamma-1.)

     iter = 0 ; temps = 0.
     do while ( temps < dt_tot)


        if (TT .lt.0) then
         write(*,*) 'prob Temp',TT, NN
!         write(*,*) 'repair assuming isobariticity'
         NN = max(NN,smallr)
         TT = min(4000./NN,8000.)  !2.*4000. / NN 
        endif


        TTold = TT       

        !NN is assumed to be in cc and TT in Kelvin
        call chaud_froid_2(TT,NN,ref,dRefdT)


!       write(*,*) 'check',TTold, TT,NN,ref,dRefdT,iter


        if (iter == 0) then
           if (dRefDT .ne. 0.) then 
              dt = abs(1.0E-1 * alpha/dRefDT) 
           else 
              dt = 1.0E-1 * dt_tot 
           endif
           dt_max = dt_tot - temps
           if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        endif

        dTemp = ref/(alpha/dt - dRefdT) 

       eps = abs(dTemp/TT)
        if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

        TT = TTold + dTemp
        if (TT < 0.) then
           write(*,*) 'Temperature negative !!!'
           write(*,*) 'TTold,TT   = ',TTold,TT
           write(*,*) 'rho   = ',rho 
           TT = 100.  !*kelvin
        endif


        iter = iter + 1

        temps = temps + dt

        dt = vardt*varrel*dt/Max(vardt*eps, varrel)

        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
!        write(*,987) temps, TT
!987     format(E10.3,2x,E10.3)
!        read (*,*)
     enddo


!  if (TT .ge. 50.)  TT=50.

!!now convert temperature in code units
TT = TT / scale_T2


return
end subroutine calc_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine chaud_froid_2(T,n,ref,dRefDT)

  use amr_parameters

  implicit none

  real(dp) :: T,P,N,x,ne                         !x est le taux d'ionisation
  real(dp) :: T2, ref2,dRefDT
  real(dp) :: froid,chaud,ref,nc, froidc    
  real(dp) :: froidcII, froido, froidh
  real(dp) :: froidc_m,froidcII_m,froido_m
  real(dp) :: param, G0, epsilon,k1,k2,bet,froidrec
  real(dp) :: eps
  

  !fonction de chauffage et de refroidissement calculee a partir des
  !refroidissements des differents elements



  !abondance de carbone 3.5 10-4, depletion 0.4

!!! on calcule l'ionisation 
!!! on suppose que si x est superieure a 1.d-4 elle est domine par 
!!! l'hydrogene et que la valeur inf est donne par le carbone 
!!! et vaut 3.5 1.d-4 * depletion * densite

!!! Pour les electrons dus a l'hydrogene on prend la 
!!! formule donnee par Wolfire et al. 2003 appendice C2.
!!! on prend un taux d'ionisation de 1.d-16 G0'=GO/1.7
!!! Z'd = 1 et phipah=0.5



  ne = 2.4E-3*((T/100.)**0.25)/0.5 !formule C15 Wolfire et al. 2003 

  x = ne / N   ! ionisation

  x = min(x,0.1)

  x = max(x,3.5E-4*0.4)


  !transition hyperfine a basse temperature: carbone et oxygene
  !chiffre pris dans la these de Karl Joulain 

  !refroidissement par le carbone




  !      froidcII = ( 2.2d-23                     !excitation par H
  !     c          + 5.5d-20 * 2 / sqrt(T) * x ) !excitation par e 
  !     c              * 3.5d-4 * 0.4d0 * exp(-92.d0 / T)


  froidcII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T/100.)**(-0.5))*x + 8.E-10*((T/100.)**(0.07))) &
       * 3.5E-4 * 0.4 * exp(-92./ T)
  !     c               3.d-4  * exp(-92. / T)


  !refroidissement par l'oxygene 
  !abondance 8.6 10-4 depletion 0.8

  froido = 1.E-26 * sqrt(T) * (24. * exp(-228./ T) + 7. * exp(-326./ T) )


  !      froido =  230.d0*1.38d-16 * (
  !     c            1.4d-8*x + 9.2d-11 *(T /100.d0)**(0.67) ) 
  !     c     * exp(-230.d0 / T)  

  !      froido = froido +  330.d0*1.38d-16 *( 
  !     c            1.4d-8*x + 4.3d-11 *(T /100.d0)**(0.8) ) 
  !     c     * exp(-330.d0 / T) 

  !      froido = froido +  98.d0*1.38d-16 * (
  !     c            5.d-9 *x + 1.1d-10* (T /100.d0)**(0.44) ) 
  !    c      * exp(-98.d0 / T) 


  !       froido = 2.5d-27 * (T/100)**0.4 * exp(-228.d0 / T)  


  !on tient compte de l'abondance du 
  froido = froido * 4.5E-4 


  !refroidissement par l'hydrogene 
  ! formule de Spitzer 1978
  froidh = 7.3E-19 * x * exp(-118400./ T )


  !refroidissement par les raies metastables des metaux
  !chiffre pris dans Hollenbach and McKee 1989 (ApJ 342, 306)





  !carbone une fois ionise ,1 transitions 2P 4P
  ! le poids est 1
  ! 2P->4P : 
  ! les expressions des coefficients d'excitation ont une dependance
  !en la temperature differente au dela de 10000K 
  !abondance 3.5 d-4 depletion 0.4

  !       froidcII_m = 6.2d4 * 1.38d-16 * 1.d0 *     !transition 2P->4P
  !     c ( 2.3d-8* (T/10000.)**(-0.5) * x + 1.d-12 ) *exp(-6.2d4 / T) 
  !     c    * 3.5d-4 * 0.4 




  !       if ( T .le. 1.d4 ) then 
  !       froido_m = 2.3d4 * 1.38d-16 / 3.d0 *   
  !     c ( 5.1d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.3d4/T) 
  !      
  !       froido_m = froido_m + 
  !     c       4.9d4 * 1.38d-16 / 3.d0  *   
  !     c ( 2.5d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-4.9d4/T) 
  !

  !       froido_m = froido_m + 
  !     c       2.6d4 * 1.38d-16 * 1.d0  *   
  !     c ( 5.2d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.6d4/T) 

  !       else 

  !       froido_m = 2.3d4 * 1.38d-16 / 3.d0 *    
  !     c ( 5.1d-9 * (T/10000.)**(0.17) * x + 1.d-12) *exp(-2.3d4/T) 
  !      
  !       froido_m = froido_m + 
  !     c       4.9d4 * 1.38d-16 / 3.d0  *     
  !     c ( 2.5d-9 * (T/10000.)**(0.13) * x + 1.d-12) *exp(-4.9d4/T) 


  !       froido_m = froido_m + 
  !     c       2.6d4 * 1.38d-16 * 1.d0  *     
  !     c ( 5.2d-9 * (T/10000.)**(0.15) * x + 1.d-12) *exp(-2.6d4/T) 


  !       endif

  !! abondance de l'oxygene
  !       froido_m = froido_m *   4.5d-4 



!!! on somme les refroidissements
  froid = froidcII  + froidh  + froido  !+ froido_m !+  froidcII_m 


  !      froid=froid*1.d-13    !conversion en MKS


  !refroidissement par le carbone neutre. On suppose l'equilibre
  ! de la reaction C + hv <=> C+ + e-
  ! les taux de reactions et de refroidissement sont pris dans
  !la these de Karl Joulain. 

  ! abondance du carbone relative a n (MKS)     


  !    C+ + e- => C 
  !       k1 = 4.4d-12 * (T/300.)**(-0.61) !s^-1 cm^-3

  !       k1 = k1 


  !    C => C+ + e-
  !       k2 = 2.2d-10 


  ! on a : [C] = k1/k2 [C+] * [e-]
  ! on suppose que tout le carbone est sous forme C+
  ! et que [e-] = [C+]

  ! l'abondance relative de carbone
  !      nc = k1/k2 * (3.5d-4*0.4)**2 * n


  !      froidc =  1.0d-24 * ( 1.4d0 * exp( -23.d0 / T )  
  !     c                + 3.8d0 * exp( -62.d0 / T )   )

  !      froidc = froidc * nc !(nc est l'abondance relative du carbone)


  !       n=exp(log(10.d0)*logn) !ici pas besoin de log

  !       valeur utilisees au 23/08/98
  !       chaud=4.d0*exp(-24.5d0*log(10.d0))*1.d-7  !un peu empirique ....



!!!! on calcule le chauffage
!!! on prend en compte le chauffage sur les grains 
!!! formules 1 et 2  de Wolfire et al. 1995 

!!!! G0 est le flux UV par rapport au flux defini par Habing et 
!!!! Draine

  G0 = 1./1.7

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T/1.E4)**0.7 / (1. + (param/5.E3) )


  chaud = 1.E-24 * epsilon 


  ! pour un flux de rayonnement G0/1.7 
  chaud = chaud * G0  

  !refroidissement recombinaison sur les grains charges positivement
  bet = 0.74/(T**0.068)
  froidrec = 4.65E-30*(T**0.94)*(param**bet)*x



  !! chaud=1.d-32 !un peu empirique ....

  !      froidc=0.d0


  ref= chaud*n - (n**2)*(froid + froidrec) !!!+ froidc)

!------------------------------------------------------

  eps = 1.0E-5
  T2 = T*(1.+eps)

  ne = 2.4E-3*((T2/100.)**0.25)/0.5 !formule C15 Wolfire et al. 2003 
  x = ne / N                       ! ionisation
  x = min(x,0.1)
  x = max(x,3.5E-4*0.4)


  froidcII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T2/100.)**(-0.5))*x + 8.E-10*((T2/100.)**(0.07))) &
       * 3.5E-4 * 0.4 * exp(-92./ T2)

  froido = 1.E-26 * sqrt(T2) * (24. * exp(-228./ T2) + 7. * exp(-326./ T2) )
  froido = froido * 4.5E-4 

  froidh = 7.3E-19 * x * exp(-118400./ T2 )

  froid = froidcII  + froidh  + froido  !+ froido_m !+  froidcII_m 

  G0 = 1./1.7

  param = G0 * sqrt(T2)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T2/1.E4)**0.7 / (1. + (param/5.E3) )


  chaud = 1.E-24 * epsilon 

  chaud = chaud * G0  

  bet = 0.74/(T2**0.068)
  froidrec = 4.65E-30*(T2**0.94)*(param**bet)*x

  ref2= chaud*n - (n**2)*(froid + froidrec) !!!+ froidc)

  dRefDT = (ref2-ref)/T/eps

  return 

end subroutine chaud_froid_2

