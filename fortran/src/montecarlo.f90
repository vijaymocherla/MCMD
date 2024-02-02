module montecarlo
    
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    implicit none

    private 

    public :: calculate_energy, update_energy
    public :: mcmove, metropolis

    !public data
    integer, public  :: n_part             ! no. of particles 
    real(dp), public :: box(3)             ! box size  
    real(dp), public :: r_cut              ! LJ cut-off distance
    real(dp), public :: e_cut              ! cut-off energy 
    real(dp), public :: temperature        ! inverse temperature
    real(dp), public :: delta_max          ! maximum allowed displacement 
    real(dp), allocatable, dimension(:,:), public :: x ! positions
    
    contains

    subroutine calculate_energy(en)
        use, intrinsic :: iso_fortran_env, only:dp=>real64  
    
        implicit none
        real(dp), intent(out) :: en
        real(dp) :: r2, r2i, r6i, xr(3)
        integer :: i, j
        en = 0.0d0
        do i = 1, n_part-1          
            do j = i+1, n_part
                xr = x(i,:) - x(j,:)        ! separation b/w (i,j)
                xr = xr - nint(xr/box)*box  ! periodic boundary conditions
                r2 = sum(xr**2)             ! sqr'd distance b/w (i,j)
                ! Check to consider interactions within cut-off distance
                if (r2 .lt. r_cut**2) then
                    r2i = 1.0d0/r2             ! inverse sqr'd distance
                    r6i = r2i**3            ! (1./r*)^6
                    en = en + 4.0d0*r6i*(r6i - 1.0) - e_cut
                end if
            end do
        end do
    end subroutine calculate_energy
    
    subroutine update_energy(i, xi_new, en)
        use, intrinsic :: iso_fortran_env, only:dp=>real64      
        implicit none
        integer, intent(in) :: i
        real(dp), intent(in) :: xi_new(3)
        real(dp), intent(inout) :: en
        real(dp) :: r2, r2i, r6i, xr(3)
        integer :: j

        do j = 1, n_part          
            if (i .ne. j) then
                ! check and remove interactions for old configuration
                xr = x(j,:) - x(i, :)            ! separation b/w (i,j)
                xr = xr - nint(xr/box)*box       ! periodic boundary conditions
                r2 = sum(xr**2)                  ! sqr'd distance b/w (i,j) 
                ! Check to consider interactions within cut-off distance
                if (r2 .lt. r_cut**2) then
                    r2i = 1.0d0/r2             ! inverse sqr'd distance
                    r6i = r2i**3            ! (1./r*)^6
                    en = en - 4.0d0*r6i*(r6i - 1.0d0) + e_cut
                end if
                ! calculate interactions in new configuration
                xr = x(j,:) - xi_new             ! separation b/w (i,j)
                xr = xr - nint(xr/box)*box       ! periodic boundary conditions
                r2 = sum(xr**2)                  ! sqr'd distance b/w (i,j)
                ! Check to consider interactions within cut-off distance
                if (r2 .lt. r_cut**2) then
                    r2i = 1.0d0/r2             ! inverse sqr'd distance
                    r6i = r2i**3            ! (1./r*)^6
                    en = en + 4.0d0*r6i*(r6i - 1.0d0) - e_cut
                end if
            end if
        end do
    end subroutine update_energy
    

    subroutine mcmove(i, decision, en)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        
        integer, intent(out)    :: i
        logical,  intent(out)   :: decision
        real(dp), intent(inout) :: en
        real(dp) :: eno, enn, delta
        real(dp) :: delx(3), xi(3), ipos
    
        ! randomly pick a particle
        n_part = size(x, dim=1)
        call random_number(ipos)
        i =  int(ipos*n_part) + 1 
        ! now randomly move the particle
        call random_number(delx)
        delx = 2.0d0*delta_max*(delx- 0.5d0)
        xi = x(i,:)
        xi = xi + delx
        xi = xi - nint(xi/box - 0.5d0)*box ! applying PBC
        ! calculate energy of the new configuration
        enn = en
        eno = en
        call update_energy(i, xi, enn)
        delta = (1.0d0/temperature)*(enn - eno) 
        call metropolis(delta, decision)
        if (decision) then
            x(i,:) = xi ! accept the move
            en = enn    ! return the new energy
        end if
    end subroutine mcmove
    
    subroutine metropolis(delta, decision)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        real(dp), intent(in) :: delta
        logical, intent(out) :: decision
    
        real(dp) :: zeta
        real(dp), parameter :: exp_cutoff = 35.0d0
    
        if (delta > exp_cutoff) then ! if too high reject
            decision = .false.
        else if (delta < 0.0d0) then   ! lowers energy accept
            decision = .true.
        else                         ! call a random number to decide
            call random_number(zeta) ! picks zeta from (0,1)
            ! metropolis test (delta = 0.0 will be accepted)  
            decision = exp(-delta) > zeta
            ! if (decision) then
            !     print*, 'GOOD LORD', delta, exp(-delta), zeta 
            ! end if
        end if
    end subroutine metropolis

end module montecarlo