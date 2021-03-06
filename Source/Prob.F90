subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module, only: probtype, p_ambient, dens_ambient, exp_energy, r_init, nsub
  use prob_params_module, only: center
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  namelist /fortin/ probtype, p_ambient, dens_ambient, exp_energy, r_init, nsub
  
  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults

  p_ambient = 1.e-5_rt        ! ambient pressure (in erg/cc)
  dens_ambient = 1.e0_rt      ! ambient density (in g/cc)
  exp_energy = 1.e0_rt        ! absolute energy of the explosion (in erg)
  r_init = 0.05e0_rt          ! initial radius of the explosion (in cm)
  nsub = 4

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set local variable defaults
  center(1) = (problo(1)+probhi(1))/2.e0_rt
  center(2) = (problo(2)+probhi(2))/2.e0_rt
  center(3) = (problo(3)+probhi(3))/2.e0_rt
  
end subroutine amrex_probinit



subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)

  use probdata_module, only: probtype, r_init, exp_energy, nsub, p_ambient, dens_ambient
  use actual_eos_module, only: gamma_const
  use bl_constants_module, only: M_PI, FOUR3RD
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only: center
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

  real(rt) :: xmin, ymin, zmin
  real(rt) :: xx, yy, zz
  real(rt) :: dist
  real(rt) :: eint, p_zone
  real(rt) :: vctr, p_exp

  integer :: i,j,k, ii, jj, kk
  integer :: npert, nambient
  
  if (probtype .eq. 32) then

     ! set explosion pressure -- we will convert the point-explosion
     ! energy into a corresponding pressure distributed throughout the
     ! perturbed volume
     vctr  = M_PI*r_init**2
     p_exp = (gamma_const - 1.e0_rt)*exp_energy/vctr
     
     do k = lo(3), hi(3)
        zmin = xlo(3) + delta(3)*dble(k-lo(3)) 
        
        do j = lo(2), hi(2)
           ymin = xlo(2) + delta(2)*dble(j-lo(2))
           
           do i = lo(1), hi(1)
              xmin = xlo(1) + delta(1)*dble(i-lo(1))
              
              npert = 0
              nambient = 0
              
              do jj = 0, nsub-1
                 yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5e0_rt)
                 
                 do ii = 0, nsub-1
                    xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5e0_rt)
                    
                    dist = (center(1)-xx)**2 + (center(2)-yy)**2
                    
                    if(dist <= r_init**2) then
                       npert = npert + 1
                    else
                       nambient = nambient + 1
                    endif
                    
                 enddo
              enddo
              
              p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient) / &
                       (dble(npert) + dble(nambient))

              eint = p_zone/(gamma_const - 1.e0_rt)

              state(i,j,k,URHO) = dens_ambient
              state(i,j,k,UMX) = 0.e0_rt
              state(i,j,k,UMY) = 0.e0_rt
              state(i,j,k,UMZ) = 0.e0_rt
              
              state(i,j,k,UEDEN) = eint +  &
                   0.5e0_rt*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMZ)**2/state(i,j,k,URHO))

              state(i,j,k,UEINT) = eint

              state(i,j,k,UFS) = state(i,j,k,URHO)

           enddo
        enddo
     enddo
     
  else if (probtype .eq. 33) then

     ! set explosion pressure -- we will convert the point-explosion energy into
     ! a corresponding pressure distributed throughout the perturbed volume
     vctr  = FOUR3RD*M_PI*r_init**3
     p_exp = (gamma_const - 1.e0_rt)*exp_energy/vctr

     do k = lo(3), hi(3)
        zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

        do j = lo(2), hi(2)
           ymin = xlo(2) + delta(2)*dble(j-lo(2))
           
           do i = lo(1), hi(1)
              xmin = xlo(1) + delta(1)*dble(i-lo(1))

              npert = 0
              nambient = 0

              do kk = 0, nsub-1
                 zz = zmin + (delta(3)/dble(nsub))*(kk + 0.5e0_rt)
                 
                 do jj = 0, nsub-1
                    yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5e0_rt)
                    
                    do ii = 0, nsub-1
                       xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5e0_rt)
                       
                       dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2
                       
                       if(dist <= r_init**2) then
                          npert = npert + 1
                       else
                          nambient = nambient + 1
                       endif
                       
                    enddo
                 enddo
              enddo
              
              p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient)/  &
                   dble(nsub*nsub*nsub)

              eint = p_zone/(gamma_const - 1.e0_rt)

              state(i,j,k,URHO) = dens_ambient
              state(i,j,k,UMX) = 0.e0_rt
              state(i,j,k,UMY) = 0.e0_rt
              state(i,j,k,UMZ) = 0.e0_rt
              
              state(i,j,k,UEDEN) = eint + &
                   0.5e0_rt*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMZ)**2/state(i,j,k,URHO))

              state(i,j,k,UEINT) = eint

              state(i,j,k,UFS) = state(i,j,k,URHO)

           enddo
        enddo
     enddo
     
  else 

     call bl_error('Dont know this probtype in initdata')

  end if
  
end subroutine ca_initdata

