! advection routines in support of method of lines integration

module mol_module

  implicit none

contains
  
#ifdef CUDA
  attributes(device) &
#endif
  subroutine mol_single_stage(time, &
                              lo, hi, domlo, domhi, &
                              dx, dt, u, h, &
                              courno, verbose)

    use advection_util_module, only: compute_cfl, divu, normalize_species_fluxes, &
                                     ht, ut, ft
    use bl_constants_module, only: ZERO, HALF, ONE, FOURTH
    use flatten_module, only: uflaten
    use riemann_module, only: cmpflx
    use ppm_module, only: ppm_reconstruct
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                   UTEMP, UEINT, UMX, GDU, GDV, GDW, &
                                   QU, QV, QW, QPRES, NQAUX

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose
    integer,  intent(in   ) :: domlo(3), domhi(3)

    real(rt), intent(in   ) :: dx(3), dt, time
    real(rt), intent(inout) :: courno

    type(ht), intent(inout) :: h
    type(ut), intent(inout) :: u

    type(ft) :: f

    integer :: ngf

    integer :: edge_lo(3), edge_hi(3)
    integer :: g_lo(3), g_hi(3)

    real(rt), parameter :: difmag = 0.1d0

    real(rt) :: div1
    integer :: i, j, k, n
    integer :: kc, km, kt, k3d

    ngf = 1

    g_lo = lo - ngf
    g_hi = hi + ngf

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             h%shk(i,j,k) = ZERO
          enddo
       enddo
    enddo

    ! Check if we have violated the CFL criterion.
    call compute_cfl(u, lo, hi, dt, dx, courno)

    ! Compute flattening coefficient for slope calculations.

    call uflaten(g_lo, g_hi, u, h)

    do n = 1, NQ

       call ppm_reconstruct(u, h, lo, hi, dx, n)

       ! Construct the interface states -- this is essentially just a
       ! reshuffling of interface states from zone-center indexing to
       ! edge-centered indexing
       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                ! x-edges

                ! left state at i-1/2 interface
                h%qm(i,j,k,n,1) = h%sxp(i-1,j,k,n)

                ! right state at i-1/2 interface
                h%qp(i,j,k,n,1) = h%sxm(i,j,k,n)

                ! y-edges

                ! left state at j-1/2 interface
                h%qm(i,j,k,n,2) = h%syp(i,j-1,k,n)

                ! right state at j-1/2 interface
                h%qp(i,j,k,n,2) = h%sym(i,j,k,n)

                ! z-edges

                ! left state at k3d-1/2 interface
                h%qm(i,j,k,n,3) = h%szp(i,j,k-1,n)

                ! right state at k3d-1/2 interface
                h%qp(i,j,k,n,3) = h%szm(i,j,k,n)

             end do
          end do
       end do

    end do

    ! Compute F^x at kc (k3d)
    f%flux => u%flux1
    call cmpflx(u, f, h, 1, lo, [hi(1)+1, hi(2), hi(3)], domlo, domhi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             h%q1(i,j,k,:) = h%qint(i,j,k,:)
          end do
       end do
    end do

    ! Compute F^y at kc (k3d)
    f%flux => u%flux2
    call cmpflx(u, f, h, 2, lo, [hi(1), hi(2)+1, hi(3)], domlo, domhi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             h%q2(i,j,k,:) = h%qint(i,j,k,:)
          end do
       end do
    end do

    ! Compute F^z at kc (k3d)

    f%flux => u%flux3
    call cmpflx(u, f, h, 3, lo, [hi(1), hi(2), hi(3)+1], domlo, domhi)

    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             h%q3(i,j,k,:) = h%qint(i,j,k,:)
          end do
       end do
    end do

    ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
    edge_lo = lo
    edge_hi = hi + 1
    call divu(lo, hi, u, h, dx)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             h%pdivu(i,j,k) = &
                  HALF*(h%q1(i+1,j,k,GDPRES) + h%q1(i,j,k,GDPRES)) * &
                       (h%q1(i+1,j,k,GDU) - h%q1(i,j,k,GDU))/dx(1) + &
                  HALF*(h%q2(i,j+1,k,GDPRES) + h%q2(i,j,k,GDPRES)) * &
                       (h%q2(i,j+1,k,GDV) - h%q2(i,j,k,GDV))/dx(2) + &
                  HALF*(h%q3(i,j,k+1,GDPRES) + h%q3(i,j,k,GDPRES)) * &
                       (h%q3(i,j,k+1,GDW) - h%q3(i,j,k,GDW))/dx(3)
          enddo
       enddo
    enddo

    do n = 1, NVAR

       if ( n == UTEMP ) then
          u%flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
          u%flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
          u%flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO

       else
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                   div1 = FOURTH*(h%div(i,j,k) + h%div(i,j+1,k) + &
                                  h%div(i,j,k+1) + h%div(i,j+1,k+1))
                   div1 = difmag*min(ZERO,div1)

                   u%flux1(i,j,k,n) = u%flux1(i,j,k,n) + &
                        dx(1) * div1 * (u%uin(i,j,k,n)-u%uin(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                   div1 = FOURTH*(h%div(i,j,k) + h%div(i+1,j,k) + &
                                  h%div(i,j,k+1) + h%div(i+1,j,k+1))
                   div1 = difmag*min(ZERO,div1)

                   u%flux2(i,j,k,n) = u%flux2(i,j,k,n) + &
                        dx(2) * div1 * (u%uin(i,j,k,n)-u%uin(i,j-1,k,n))
                enddo
             enddo
          enddo

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   div1 = FOURTH*(h%div(i,j,k) + h%div(i+1,j,k) + &
                                  h%div(i,j+1,k) + h%div(i+1,j+1,k))
                   div1 = difmag*min(ZERO,div1)

                   u%flux3(i,j,k,n) = u%flux3(i,j,k,n) + &
                        dx(3) * div1 * (u%uin(i,j,k,n)-u%uin(i,j,k-1,n))
                enddo
             enddo
          enddo

       endif

    enddo

    call normalize_species_fluxes(u, lo, hi)

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                u%update(i,j,k,n) = u%update(i,j,k,n) + &
                     (u%flux1(i,j,k,n) * u%area1(i,j,k) - u%flux1(i+1,j,k,n) * u%area1(i+1,j,k) + &
                      u%flux2(i,j,k,n) * u%area2(i,j,k) - u%flux2(i,j+1,k,n) * u%area2(i,j+1,k) + &
                      u%flux3(i,j,k,n) * u%area3(i,j,k) - u%flux3(i,j,k+1,n) * u%area3(i,j,k+1) ) / u%vol(i,j,k)

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   u%update(i,j,k,n) = u%update(i,j,k,n) - h%pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                u%flux1(i,j,k,n) = dt * u%flux1(i,j,k,n) * u%area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                u%flux2(i,j,k,n) = dt * u%flux2(i,j,k,n) * u%area2(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                u%flux3(i,j,k,n) = dt * u%flux3(i,j,k,n) * u%area3(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine mol_single_stage

end module mol_module
