! diffusion.f90
!
! Adapted from  Allen, M. P., & Tildesley, D. J. (2017). Computer simulation of liquids. Oxford university press. 
! Source: https://github.com/Allen-Tildesley/examples
! 
program diffusion
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use io_module

    implicit none
    interface
    subroutine unfold(x_old, x, box)
        double precision, dimension(:,:), intent(in)    :: x_old ! Previous configuration
        double precision, dimension(:,:), intent(inout) :: x     ! Current configuration
        double precision, intent(in) :: box(3)    
    end subroutine
    end interface
    integer :: n_part          ! Number of atoms
    integer :: nt              ! Number of timesteps to correlate
    integer :: n0              ! Number of time origins to store
    integer :: dt0             ! Interval for time origins
    integer :: dt              ! Time difference (in timesteps)
    integer :: t               ! Time (in timesteps, equivalent to number of file)
    real(dp) :: delta          ! Time interval (simulation units: only used in output file)
  
    real(dp) :: sigma_angs
    real(dp) :: angs_per_ps
    real(dp) :: time_ps
    real(dp) :: box(3)
    
    real(dp), dimension(:,:),   allocatable :: x     ! Positions (n_part, 3)
    real(dp), dimension(:,:),   allocatable :: x_old ! Previous positions (n_part, 3)
    real(dp), dimension(:,:),   allocatable :: v     ! Velocities (n_part, 3)
    real(dp), dimension(:,:,:), allocatable :: x0    ! Stored position origins (n_part, 3, n0)
    real(dp), dimension(:,:,:), allocatable :: v0    ! Stored velocity origins (n_part, 3, n0)
    real(dp), dimension(:),     allocatable :: vacf  ! Velocity correlation function (0:nt)
    real(dp), dimension(:),     allocatable :: norm  ! array to keep track of normalisation
    real(dp), dimension(:),     allocatable :: xvcf  ! Displacement-velocity cross correlation function (0:nt)
    real(dp), dimension(:),     allocatable :: msd   ! Mean-squared displacement (0:nt)
    integer,  dimension(:),     allocatable :: t0    ! Times of origins

    integer :: k, mk, nk
    integer :: nstart
    integer :: rc, ci
    logical :: full
    character (len=1024)  :: xyz_file
    character (len=1024)  :: velocity_file
    character (len=1024)  :: outfile
    character (len=1024)  :: comment
    character (len=1024)  :: arg
    character (len=32)    :: boxvar
    character (len=2)     :: element

    ! Example default values
    nt              = 50   ! Max correlation time (as a multiple of interval between configurations)
    dt0             = 10   ! We only take time origins at these intervals, for efficiency
    delta           = 0.1  ! This should be set to the actual time interval between configurations
    sigma_angs      = 3.405d0
    time_ps         = 2.18d0
    angs_per_ps     = sigma_angs / time_ps 
    ci = 1  
    nstart = 0

    do
        call get_command_argument(ci, arg) 
        if (trim(arg) == "-xyz") then               ! xyz filename
            call get_command_argument(ci+1, arg) 
            xyz_file = trim(arg)
            ci = ci + 2
        else if (trim(arg) == "-vel") then          ! velocity filename
            call get_command_argument(ci+1, arg) 
            velocity_file = trim(arg)
            ci = ci + 2    
        else if (trim(arg) == "-o") then            ! output filename
            call get_command_argument(ci+1, arg) 
            outfile = trim(arg)
            ci = ci + 2
        else if (trim(arg) == "-nt") then           ! max correlation time    
            call get_command_argument(ci+1, arg) 
            read(arg,'(i16)') nt
            ci = ci + 2  
        else if (trim(arg) == "-dt0") then          ! origin_interval    
            call get_command_argument(ci+1, arg) 
            read(arg,'(i16)') dt0
            ci = ci + 2  
        else if (trim(arg) == "-delta") then        ! input conversion factor    
            call get_command_argument(ci+1, arg) 
            read(arg,'(f32.16)') delta
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

    
    open(100, file=xyz_file)
    open(200, file=velocity_file)
    call read_xyz(x, sigma_angs, comment, element, 100, rc)
    if (is_iostat_end(rc)) then
        print*, "Error in reading velocity file"
        stop
    end if
    call read_xyz(v, angs_per_ps, comment, element, 200, rc)
    if (is_iostat_end(rc)) then
        print*, "Error in reading velocity file"
        stop
    end if

    ! initialising variables
    n_part = size(v, dim=1)
    read(comment,*)  boxvar, box
    box = box/sigma_angs
    n0 = nt / dt0 + 1 ! Enough origins to span max correlation time
    allocate ( x_old(n_part, 3), x0(n_part,3,n0), v0(n_part,3,n0), t0(n0) )
    allocate ( msd(0:nt), xvcf(0:nt), vacf(0:nt), norm(0:nt) )
    
    msd  = 0.0d0
    xvcf = 0.0d0
    vacf = 0.0d0
    t    = 0 ! Time of each snapshot of the trajectory
    mk   = 0 ! Storage location of time origin
    full = .false.
    ! storing time-step at t=0
    mk = mk + 1
    x0(:,:,mk) = x(:,:) ! Store position origins
    v0(:,:,mk) = v(:,:) ! Store velocity origins

    do ! Loop reading and correlating data
        call read_xyz(x, sigma_angs, comment, element, 100, rc)
        call read_xyz(v, angs_per_ps, comment, element, 200, rc)
        if (is_iostat_end(rc)) exit ! terminate when EoF is reached
        
        ! remove the effect of PBC
        call unfold ( x_old, x, box )
        
        if ( modulo( t, dt0 ) == 0 ) then ! Test to store as time origin
                mk = mk + 1 ! Increment origin counter
    
                if ( mk > n0 ) then ! Test for over-running origin arrays
                full = .TRUE.
                mk = mk - n0 ! Wrap-around to overwrite older origins
            end if ! end test for over-running origin arrays
            
                t0(mk)     = t      ! Store time origin
                x0(:,:,mk) = x(:,:) ! Store position origins
                v0(:,:,mk) = v(:,:) ! Store velocity origins
                
            end if ! end test to store as time origin
            
            if ( full ) then
                nk = n0 ! Correlate with all stored time origins
            else
                nk = mk ! Correlate with those stored so far
            end if
    
            do k = 1, nk ! Loop over time origins
                
                dt = t - t0(k)  ! Time difference between current configuration and time origin
                
                
                if ( dt <= nt ) then ! Check that dt is within desired range
                    
                    ! Sum over all N atoms and over 3 Cartesian components
                    msd(dt)  = msd(dt)  + sum(( x(:,:) - x0(:,:,k) ) ** 2)        ! Increment msd
                    xvcf(dt) = xvcf(dt) + sum(( x(:,:) - x0(:,:,k) ) * v0(:,:,k)) ! Increment cross correlation function
                    vacf(dt) = vacf(dt) + sum(  v(:,:) * v0(:,:,k))               ! Increment autocorrelation function    
                    norm(dt) = norm(dt) + 1.0
                end if ! end check that dt is within desired range
                
            end do ! end loop over time origins    
            
            x_old = x     ! Ready to unfold next step
            t = t + 1 ! Number of next step

    end do ! end loop reading and correlating data
    close(100)
    close(200)

    ! Normalize by N as well as time-origin normalizing factors
    msd  = msd  / norm / n_part ! 3D mean-squared displacement
    xvcf = xvcf / norm / n_part ! 3D cross-correlation function
    vacf = vacf / norm / n_part ! 3D autocorrelation function

    delta = delta*time_ps
    open(10,file=outfile)    
    write(10,'(i16)') nt 
    write(10,'(1a16,3a32)') 'Tau (ps)', '<v(t) . v(t+T)>', '<(x(t)-x(t+T)) . v(t+T)>', '<(x(t)-x(t+T))**2>'
    do t = 0, nt
        write (10, fmt='(f16.8,3f32.16)' ) t*delta, vacf(t), xvcf(t), msd(t) ! Include delta in time
    end do
    close(10)

end program diffusion

subroutine unfold (x_old, x, box)
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    implicit none
    real(dp), dimension(:,:), intent(in)    :: x_old ! Previous configuration
    real(dp), dimension(:,:), intent(inout) :: x     ! Current configuration
    real(dp), intent(in) :: box(3)
    integer :: n, i
    ! Removes effects of periodic boundaries on particle trajectories
    ! Assumes that no particle actually moves more than half a box length
    ! between successive configurations

    n = size(x, dim=1)
    x = x - x_old                   ! Convert r to displacements relative to r_old
    x = x - x_old                   ! Convert r to displacements relative to r_old
    do i=1,n                        ! Apply periodic boundaries to displacements
        x(i,:) = x(i,:) - nint( x(i,:) / box ) * box   
    end do
    x = x + x_old                   ! Convert r back to absolute coordinates

end subroutine unfold
