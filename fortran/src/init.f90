! initialize.f90
! Sets up an initial configurations and velocites for MC/MD simulations.
!
!   Written by Sai Vijay Mocherla <vijaysai.mocherla@gmail.com>


program main
    ! Generates initial configuration of atoms (positions, velocities)
    ! Reads input from init.inp in standard namelist &init_params 
    
    ! Things to do:
    !  - Make sure everything is in reduced units                          
    !  - Setup velocity sampling from boltzmann distribution
    !  - Calculate melting factor (S0) to check lattices.                
    
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use initialize, only : get_positions, get_velocities, check_overlap, &
                                  allocate_arrays, deallocate_arrays, &
                                  n_part, x, v  
    use io_module, only : write_xyz

    implicit none
    
    integer              :: rc            ! iostat variable 
    real(dp)             :: rho           ! particle number density 
    real(dp)             :: mass          ! mass of the particles
    real(dp)             :: temperature   ! ensemble temperatureerature
    real(dp)             :: box(3)        ! size of the simulation box
    real(dp)             :: sigma_angs    ! sigma in angs for conversion
    real(dp)             :: angs_per_ps   ! velocity conversion factor  
    real(dp)             :: time_ps       ! time conversion factor
    ! string variable
    character (len=2)    :: element       ! atom type   
    character (len=1024) :: xyz_file      ! file to write positions
    character (len=1024) :: velocity_file ! file to write velocities
    character (len=1024) :: infile        ! initialize input file (w/ std. nml)   
    character (len=1024) :: outfile       ! initialize output file   

    character (len=1024) :: sample_type   ! v-sampling: 'uniform', 'boltzmann'
    character (len=1024) :: fill_type     ! fills on a lattice or randomly
    character (len=1024) :: comment       ! comment line in .xyz and .vxyz files
    logical              :: positions     ! switch for generating positions
    logical              :: velocities    ! switch for generating velocities 
    ! cli variables
    integer              :: ci
    character (len=1024) :: arg
    namelist /init_params/ positions, velocities, element, n_part, &
                           xyz_file, rho, sigma_angs, fill_type, &
                           velocity_file, mass, temperature, sample_type
    ! Set default parameters
    positions     = .true.                 ! positions are always generated
    velocities    = .false.                ! initial velocities are generated when needed 
    element       = 'Ar'                   ! by default, we use argon atoms 
    n_part        = 216                    ! we consider (nc  =6)**3 particles
    xyz_file      = 'cubic_216.xyz'        ! positions are stored in cubic_216.xyz
    rho           = 0.60d0                 ! should lie in the liquid region for argon at T*>0.728
    fill_type     = 'cubic'                ! simple cubic is always considered
    velocity_file = 'boltzmann_216.vxyz'   ! positions are stored in cubic_216.xyz
    mass          = 1.0d0                  ! we scale all masses with that of an argon atom
    temperature   = 1.0d0                  ! should lie in the liquid region for argon
    sample_type   = 'boltzmann'            ! samples from a boltzmann distribution
    sigma_angs    = 3.405d0

    ! Getting command line input
    ci = 1
    do
        call get_command_argument(ci, arg) 
        if (trim(arg) == "-i") then
            call get_command_argument(ci+1, arg) 
            infile = trim(arg)
            ci = ci + 2
        else if (trim(arg) == "-o") then
            call get_command_argument(ci+1, arg) 
            outfile = trim(arg)
            ci = ci + 2    
        else
            exit
        end if
    end do

    ! Get input parameters from init.inp
    open(unit=10, file=infile)
    read(unit=10, nml=init_params, iostat=rc)
    if (rc /= 0) then
        print*, '!! ERROR: Can not read input file.'
        stop
    end if
    close(10)
    
    ! calculating conversion velocity factor
    time_ps       = 2.180d0
    angs_per_ps   = sigma_angs/time_ps
    ! Calculating dimensions of the box assuming it is cube. 
    box(1) = (n_part/rho)**(1.0d0/3.0d0)
    box(2) = box(1)
    box(3) = box(1)
    open(unit=200, file=outfile)
    write(200, fmt='(a,a)') 'Initialize'
    write(200, *)
    
    ! Generate initial positions
    call allocate_arrays(positions, velocities)
    if (positions) then
        write(200, *)
        write(200, *) 'Generating Initial Positions:'
        write(200, *) 'Box Size (lx, ly, lz) in Angstrom = ', sigma_angs*box
        write(comment, fmt='(a32,3f16.8)') "'!  box_size (in angstrom)  =  '", sigma_angs*box(:)
        write(200, *) 'Placing particles in ',  trim(fill_type), ' positions.'
        call get_positions(box, fill_type)
        write(200, *) 'Writing co-ordinates to ', trim(xyz_file)
        open(unit=300, file=xyz_file)
            call write_xyz(x, sigma_angs, comment, element, unit=300)
        close(300)
        ! checking if there are any overlapping particles
        call check_overlap
        write(200, *) 'Done!'
        write(200, *)
    end if
    ! Generate initial velocity distribution
    if (velocities) then
        write(200, *)
        write(200, *) 'Generating Initial Velocities:'
        write(200, *) 'Box Size (lx, ly, lz) = ', sigma_angs*box
        write(comment, fmt='(a32,3f16.8)') "'!  box_size (in angstrom)  =  '", sigma_angs*box(:)
        write(200, *) 'Sampling velocities from ',  trim(sample_type), ' distribution.'
        call get_velocities(mass, temperature, sample_type)
        write(200, *) 'Writing co-ordinates to ', trim(velocity_file)
        write(200,'(a,f8.4)') 'Temperature of the ensemble', mass*sum(v**2)/(3*n_part)
        open(300, file=velocity_file)
            call write_xyz(v, angs_per_ps, comment, element, unit=300)
        close(300)
        write(200, *) 'Done!'
        write(200, *)
    end if
    call deallocate_arrays(positions,velocities)
    close(200)

end program main

