module initialize

    use, intrinsic :: iso_fortran_env, only:dp=>real64
    implicit none

    private

    public :: allocate_arrays, deallocate_arrays
    public :: get_positions, get_velocities, check_overlap
    !public data
    integer, public :: n_part                          ! no. of particles 
    real(dp), allocatable, dimension(:,:), public :: x ! positions
    real(dp), allocatable, dimension(:,:), public :: v ! velocities
    real(dp), parameter :: pi_dp = 4.0d0*atan(1.0d0)   ! PI upto double precision value 
    real(dp), parameter :: sigma = 1.0d0
    contains

    subroutine allocate_arrays(positions, velocities)
        implicit none
        logical, intent(in) :: positions
        logical, intent(in) :: velocities
        
        if (positions) then
            allocate(x(n_part, 3))
        end if

        if (velocities) then
            allocate(v(n_part, 3))
        end if
    end subroutine

    subroutine deallocate_arrays(positions, velocities)

        implicit none
        logical, intent(in) :: positions
        logical, intent(in) :: velocities
                
        if (positions) then
            deallocate(x)
        end if

        if (velocities) then
            deallocate(v)
        end if
    end subroutine

    subroutine get_positions(box, fill_type)
        ! Returns an array of positions based on given input method
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        real(dp), intent(in) :: box(3)
        character (len=1024), intent(in) :: fill_type

        real(dp), allocatable, dimension(:,:) :: rvecs
        integer :: i
        integer :: nx, ny, nz, cells(3)
        integer :: nseed

        if (fill_type == 'random') then
            ! random fill
            nseed = 1000
            call random_seed(size=nseed)
        
            i = 1
            do while(i<=n_part)
                call random_fill(i, box)
                i = i + 1
            end do
    
        else if (fill_type == 'cubic') then
            ! 'cubic'
            nx = nint((n_part)**(1.0d0/3.0d0))
            ny = nx
            nz = nx
            cells = (/nx, ny, nz/)
            allocate(rvecs(1,3))
            rvecs(1, :)  = (/0.25d0, 0.25d0, 0.25d0/)
            call lattice_fill(cells, box, rvecs)
    
        else if (fill_type == 'bcc') then
            ! 'body-centered cubic'
            nx = nint((n_part/2.0d0)**(1.0d0/3.0d0))
            ny = nx
            nz = nx
            cells = (/nx, ny, nz/)
            allocate(rvecs(2,3))
            rvecs(1, :)  = (/0.25d0, 0.25d0, 0.25d0/)
            rvecs(2, :)  = (/0.75d0, 0.75d0, 0.75d0/)
            call lattice_fill(cells, box, rvecs)
    
        else if (fill_type == 'fcc') then
            !'face centered cubic'
            nx = nint((n_part/4.0d0)**(1.0d0/3.0d0))
            ny = nx
            nz = nx
            cells = (/nx, ny, nz/)
            allocate(rvecs(4,3))
            rvecs(1, :)  = (/0.25d0, 0.25d0, 0.25d0/)
            rvecs(2, :)  = (/0.25d0, 0.75d0, 0.75d0/)
            rvecs(3, :)  = (/0.75d0, 0.25d0, 0.75d0/)
            rvecs(4, :)  = (/0.75d0, 0.75d0, 0.25d0/)
            call lattice_fill(cells, box, rvecs)
    
        else
            print*, "!!ERROR : Unknown initialization method"
            stop
        end if
    end subroutine get_positions


    recursive subroutine random_fill(i, box)
        ! Recursively samples the box by avoiding any overlaps
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        integer, intent(in) :: i
        real(dp), intent(in) :: box(3)
        integer :: j 
        real(dp) :: distance, dr(3)
        ! random sampling from a uniform distribution
        call random_number(x(i,:))
        x(i,:) = x(i,:)*box
        ! checking for overlaps
        do j=1,i-1
            dr = x(i,:)-x(j,:)
            ! accounting for pbc
            ! PBC:
            !   o j'  |   o i -------- o j   |    o i'
            !   L = d_ij' + d_ji'
            !   if the distance > 0.5*L then 
            !       (adjust the distance for PBC)
            dr = dr - nint(dr/box)*box
            distance = dsqrt(sum((dr)**2))
            if (distance < sigma) then
                ! recursion in case of an overlap
                call random_fill(i, box)
            end if
        end do
    end subroutine random_fill


    subroutine lattice_fill(cells, box, rvecs)
        ! For a given unit cell generates the coordinates for a lattice 
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        real(dp), intent(in) :: box(3)
        integer, intent(in) :: cells(3)
        real(dp), allocatable, intent(in) :: rvecs(:,:)
        real(dp) :: r(3)
        integer :: i, j, k, l, m  ! all iter variables
        integer :: n, nx, ny, nz ! no. of cells along x,y,z
        real(dp) :: ax, ay, az ! primitive cell dimensions
        integer :: nr ! no. of position vectors for a given unit cell
        nr = size(rvecs, dim=1)
        nx = cells(1)
        ny = cells(2)
        nz = cells(3)
        ax = box(1)/cells(1)
        ay = box(2)/cells(2)
        az = box(3)/cells(3)
        m = 0 ! iterval to track no. of particles
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    do l=1,nr
                        r(1) = ax*dfloat(i-1) + ax*rvecs(l,1) 
                        r(2) = ay*dfloat(j-1) + ay*rvecs(l,2) 
                        r(3) = az*dfloat(k-1) + az*rvecs(l,3)
                        m = m + 1 
                        x(m,:) = r
                    end do
                end do
            end do
        end do
        
        if (m .ne. n_part) then
            print*, '!! ERROR : Issue with no. of particles (calculated, given): ', m , n  
            stop
        end if
    end subroutine lattice_fill


    subroutine get_velocities(mass, temperature, sample_type)
        ! Returns an array of velocities for a temperature(T) sampled using a given method.
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        real(dp), intent(in) :: mass
        real(dp), intent(in) :: temperature
        character (len=1024), intent(in) :: sample_type
        real(dp) :: kb, alpha, v_com(3)
        
        kb = 1.0 ! 1.380649E-23 ! J K^-1
        ! setting default values for optionals
        ! (not a very elegant way!)  

        call random_number(v)
        if (sample_type == 'uniform') then
            v = v - 0.5
        else if (sample_type == 'boltzmann') then
            ! sampling v in (-4*sigma, 4*sigma)
            v = 8*v 
            v = 1/(dsqrt(2.0d0*pi_dp)) * exp((v)**2 / 2)            
        end if    
        v_com = sum(v, dim=1)/n_part
        ! Correcting velocities to preserve conservation of momentum.
        v(:,1) =(v(:,1)-v_com(1))
        v(:,2) =(v(:,2)-v_com(2))
        v(:,3) =(v(:,3)-v_com(3))
        ! Scaling velocities to given temperature
        alpha = dsqrt((3*n_part*kb*temperature)/(mass*sum(v**2)))
        v = alpha*v
    end subroutine


    subroutine check_overlap
        ! Checks for overlaps b/w particles
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none        
        integer :: i, j
        real(dp) :: dr(3), distance
        do i=1,n_part-1
            do j=i+1,n_part
                dr = x(i,:) - x(j,:)
                distance = dsqrt(sum(dr**2))
                if (distance < sigma) then
                    print*, '!!ERROR: particles too close', i, j, distance
                end if
            end do
        end do
        print*, 'All checks passed!'
    end subroutine check_overlap    

end module initialize