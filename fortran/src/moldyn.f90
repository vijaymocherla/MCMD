module moldyn
    
    use, intrinsic :: iso_fortran_env, only:dp=>real64

    implicit none

    private 

    public :: force, verlet, velocity_verlet, calculate_potential
    
    ! public data
    integer,  public :: n_part             ! no. of particles 
    real(dp), public :: box(3)             ! box size  
    real(dp), public :: mass               ! mass of the particles
    real(dp), public :: r_cut              ! LJ cut-off distance
    real(dp), public :: e_cut              ! cut-off energy 
    real(dp), public :: temperature        ! inverse temperature
    real(dp), allocatable, dimension(:,:), public :: x  ! positions
    real(dp), allocatable, dimension(:,:), public :: v  ! velocities
    real(dp), allocatable, dimension(:,:), public :: f  ! forces
    real(dp), allocatable, dimension(:,:), public :: fm ! forces at ti- delt
    real(dp), allocatable, dimension(:,:), public :: xm ! forces at ti- delt
    real(dp), allocatable, dimension(:,:), public :: xx ! temporary array for positions

    contains

    subroutine force
        use, intrinsic :: iso_fortran_env,only:dp=>real64
        implicit none
      
        real(dp) :: xr(3), ff(3), r2
        integer  :: i, j    
        real(dp) :: r2i, r6i
        
        ! initialize f, en
        f = 0.0        
     
        ! loop over pairs of atoms
        do i = 1, n_part-1          
            do j = i+1, n_part
                xr = x(i,:) - x(j,:)        ! separation b/w (i,j)
                xr = xr - nint(xr/box)*box  ! periodic boundary conditions
                r2 = sum(xr**2)             ! sqr'd distance b/w (i,j)
                
                ! Check to consider interactions within cut-off distance
                if (r2 .lt. r_cut**2) then
                    r2i = 1./r2                       ! inverse sqr'd distance
                    r6i = r2i**3                      ! (1./r*)^6
                    ff = 48.0d0*r2i*r6i*(r6i - 0.5)  ! gradient of potential  
                    f(i,:) = f(i,:) + ff*xr
                    f(j,:) = f(j,:) - ff*xr                        
                end if
    
            end do
        end do
    end subroutine force

    subroutine calculate_potential(en)
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
                    r2i = 1./r2             ! inverse sqr'd distance
                    r6i = r2i**3            ! (1./r*)^6
                    en = en + 4.0d0*r6i*(r6i - 1.0) - e_cut
                end if
            end do
        end do
    end subroutine calculate_potential


    subroutine verlet(del_t) 
        use, intrinsic :: iso_fortran_env,only:dp=>real64       
        implicit none
        real(dp), intent(in) :: del_t
        integer :: i
        xx = 2*x - xm + (f/mass)*del_t**2
        v = (xx - xm)/(2*del_t)       ! update velocities to v(t + del_t)
        xm = x                        ! update prev. positions to xm(t)   
        ! Note: apply PBC after calculating velocities
        ! applying periodic boundary conditions
        do i=1,n_part
            xx(i,:) = xx(i,:) - nint(xx(i,:)/box - 0.5)*box 
        end do
        x = xx                        ! update positions to x(t + del_t)   
        call force                    ! update forces to f(t + del_t)
    end subroutine verlet


    subroutine velocity_verlet(del_t)
        use, intrinsic :: iso_fortran_env,only:dp=>real64       
        implicit none
        real(dp), intent(in) :: del_t
        integer :: i
        x = x + del_t*v + 0.5*f*del_t**2/mass  ! update positions to x(t+del_t)
        ! applying periodic boundary conditions
        do i=1,n_part
            x(i,:) = x(i,:) - nint(x(i,:)/box - 0.5)*box 
        end do
        fm = f                              ! store forces at f(t)
        call force                          ! update forces f(t + del_t) using r(t + del_t)
        v = v + 0.5*(f + fm)*del_t/mass     ! update velocities v(t + del_t) using f(t), f(t + del_t) 
    end subroutine
    

end module moldyn