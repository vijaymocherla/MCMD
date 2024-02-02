program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use io_module, only: read_xyz

    implicit none

    integer  :: n_part              ! no. of particles 
    real(dp) :: sigma_angs          ! conv. factor b/w reduced and normal units
    real(dp) :: r, xr(3)
    real(dp) :: r_lb, r_ub, vol_const
    real(dp) :: rho
    real(dp) :: box(3)              ! box size  
    real(dp) :: bw                  ! bin-width
    real(dp) :: n_r                 ! volume within a bin-width      
    integer  :: nbins               ! no. of bins
    integer  :: nstep               ! no. of time steps
    integer  :: i, j, k, ci
    integer  :: rc
    integer  :: nstart

    character (len=1024) :: xyz_file
    character (len=1024) :: outfile
    character (len=1024) :: comment
    character (len=32)   :: boxvar
    character (len=2)    :: element
    character (len=128)  :: arg
    
    real(dp), allocatable :: x(:,:) ! positions
    real(dp), allocatable :: h(:)   ! histogram of particle distances
    real(dp), allocatable :: g(:)   ! pair distribution functions
    real(dp), parameter   :: pi_dp = 4.0d0*atan(1.0d0)   ! PI upto double precision value 
    
    ci = 1
    bw = 0.02d0
    sigma_angs = 3.405d0
    nstart = 0
    ! get input arguments from commandline    
    do
        call get_command_argument(ci, arg) 
        if (trim(arg) == "-xyz") then               ! .xyz filename
            call get_command_argument(ci+1, arg) 
            xyz_file = trim(arg)
            ci = ci + 2
        else if (trim(arg) == "-o") then            ! .rdf output filename
            call get_command_argument(ci+1, arg) 
            outfile = trim(arg)
            ci = ci + 2    
        else if (trim(arg) == "-bw") then           ! input bin-width  
            call get_command_argument(ci+1, arg) 
            read(arg,'(f32.16)') bw
            ci = ci + 2  
        else if (trim(arg) == "-sigma_angs") then   ! input conversion factor    
            call get_command_argument(ci+1, arg) 
            read(arg,'(f32.16)') sigma_angs
            ci = ci + 2  
        else if (trim(arg) == "-nstart") then       ! start sampling after 'n' steps
            call get_command_argument(ci+1, arg) 
            read(arg,'(i16)') nstart
            ci = ci + 2  
        else 
            exit
        end if
    end do
    
    ! sampling trajectory file
    print*, "Sampling the trajectory file...."
    open(100,file=xyz_file)
    
    ! reading xyz at t=0
    call read_xyz(x, sigma_angs, comment, element, unit=100, rc=rc)
    
    ! initialising variables
    nstep = 0
    n_part = size(x, dim=1)
    read(comment,*)  boxvar, box
    box = box/sigma_angs
    nbins = floor(0.5d0*box(1)/bw)
    allocate(h(nbins))
    allocate(g(nbins))
    h = 0.0d0
    g = 0.0d0
 
    do  ! loop over trajectory file
        if (nstep >= nstart) then
            do i=1,n_part-1
                do j=i+1, n_part
                    xr = x(i,:) - x(j,:)           ! separation b/w (i,j)
                    xr = xr - box*nint(xr/box)     ! periodic boundary conditions
                    r = dsqrt(sum(xr**2))          ! distance b/w (i,j)
                    k = nint(r/bw) + 1
                    if (k <= nbins) then
                        h(k) = h(k) + 2   ! counting pairs
                    end if
                end do
            end do
        end if

        call read_xyz(x, sigma_angs, comment, element, unit=100, rc=rc)
        if(is_iostat_end(rc)) exit  ! Terminate when EoF is reached.
        nstep = nstep + 1
    end do 
    close(100)
    
    ! constants for calculating n_id(r)
    rho = real(n_part, 8)/box(1)**3
    vol_const = 4.0d0/3.0d0 * pi_dp * rho 

    ! calculating the pair distribution from histogram
    do k=1,nbins
        g(k) = h(k) / real(n_part*(nstep-nstart), 8)      ! averaging over timesteps and particles    
        r_lb = real(k-1, 8) * bw                 ! lower bound of the bin   
        r_ub = r_lb + bw                         ! upper bound of the bin
        n_r  = vol_const * (r_ub**3 - r_lb**3)    ! no. of particles in shell of thickness bw at 'r'
        g(k) = g(k) / n_r                        ! n(r) / n_id(r)
    end do 
    
    ! writing g(r) to outfile
    open(200, file=outfile)
    write(200, *) nbins, bw
    write(200,'(2a32)') 'r', 'g(r)'
    do k=1,nbins
        write(200,'(2f32.16)') (real(k,8)-0.5d0)*bw, g(k) 
    end do
    close(200)

    print*, 'Completed calculating pair-distribution function!'

end program main