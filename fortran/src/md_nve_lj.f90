! gfortran -O3 md.f90 initialize.f90 -o md.x
!
program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use io_module, only: read_xyz, write_xyz
    use moldyn
    implicit none
    
    ! variables for input parameters
    real(dp) :: ti              ! initial time
    real(dp) :: tf              ! final time
    real(dp) :: del_t           ! time_step (\delta t)
    integer  :: print_nstep     ! print every nsteps
    real(dp) :: sigma_angs      ! length coversion factor 
    real(dp) :: eps_kjpm        ! energy conversion factor
    real(dp) :: angs_per_ps     ! velocity conversion factor  
    real(dp) :: time_ps         ! time conversion factor

    character(len=1024) :: initial_velocities  ! initial configuration of atoms
    character(len=1024) :: initial_positions   ! initial configuration of atoms
    character(len=1024) :: outfile             ! output file to print energies  
    character(len=1024) :: xyz_file            ! xyz file to store position traj.  
    character(len=1024) :: velocity_file       ! vxyz file to store velocities traj.  
    character(len=1024) :: integration         ! integration integration
    character(len=1024) :: comment             ! comment line  
    character(len=2)    :: element             ! atom type  

    ! input from standard namelist
    namelist /md_params/ initial_positions, initial_velocities, &
                         outfile, xyz_file, velocity_file, & 
                         integration, ti, tf, del_t, print_nstep, temperature, &
                         sigma_angs, eps_kjpm, r_cut
    
    integer  :: nstep        ! iter variables
    integer  :: rc              ! iostat variables
    real(dp) :: t               ! time (iter variable)
    real(dp) :: ke              ! kinetic energy
    real(dp) :: en              ! potential energy
    real(dp) :: rc2i, rc6i
    character (len=32) :: boxvar 

    ! reading input parameter for running MD trajectory
    open(10, file='md_lj.inp')
        read(10, nml=md_params, iostat=rc)
        if (rc /= 0) then
            print*, '!! I/O ERROR: Can not read input file. IOSTAT=',rc
            stop
        end if
    close(10)
        
    ! calculating conversion factor for velocity
    time_ps = 2.180d0
    angs_per_ps = sigma_angs/time_ps
    ! initialising variables
    en = 0.0d0
    mass = 1.0d0
    ! calculating cutoff energy
    rc2i = 1.0d0 / r_cut**2
    rc6i = rc2i**3
    e_cut = 4.0d0*rc6i*(rc6i - 1.0)
    ! reading initial positions and velocities 
    open(100, file=initial_positions)
    call read_xyz(x, sigma_angs, comment, element, unit=100, rc=rc)
    read(comment, *)  boxvar, box
    box = box/sigma_angs
    close(100)
    open(100, file=initial_velocities)
    call read_xyz(v, angs_per_ps, comment, element, unit=100, rc=rc)
    close(100)
    
    ! allocating arrays
    n_part = size(x ,dim=1)
    allocate(xx(n_part,3))
    allocate(f(n_part,3))
    f = 0.0d0
    ! Time propagation 
    t = 0.0d0 
    nstep = 0

    ! opening files for i/o
    open(200,file=xyz_file)
    open(300,file=velocity_file)
    open(20, file=outfile)    
    ! writing headers for outfile
    write(20,'(a,f8.4)') 'Temperature of the ensemble (at ti)', mass*sum(v**2)/(3*n_part)
    write(20,*) int((tf-ti)/del_t)/print_nstep ! no. of time steps
    write(20,'(a16,4a32)') 'time (ps)', 'Potential Energy (KJ/mol)', 'Kinetic Energy (KJ/mol)', 'Total Energy (KJ/mol)'
    ! calculate KE and PE at t=ti, write to output
    call force  
    call calculate_potential(en)
    ke =  0.5*mass*sum(v**2)
    write(20,'(f16.8, 4f32.16)') t*time_ps, en*eps_kjpm, ke*eps_kjpm, (en+ke)*eps_kjpm
    ! writing initial coordinates to xyzfile.
    write(comment, '(a32,3f16.8)') "'!  box_size (in angstrom)  =  '", box*sigma_angs
    call write_xyz(x, sigma_angs, comment, element, unit=200)
    call write_xyz(v, angs_per_ps, comment, element, unit=300)
    
    if (integration == 'verlet') then
        allocate(xm(n_part,3))
        xm  = x - del_t*v
        do while(t <= tf)
            call verlet(del_t)
            t = t + del_t
            nstep = nstep + 1            
            if ( modulo(nstep,print_nstep) .eq. 0 ) then
                call calculate_potential(en)
                write(20,'(f16.8, 4f32.16)') t*time_ps, en*eps_kjpm, ke*eps_kjpm, (en+ke)*eps_kjpm
                ! writing initial coordinates to xyzfile.
                write(comment, '(a,f16.8)') "Time t (ps) = ", t*time_ps
                call write_xyz(x, sigma_angs, comment, element, unit=200)
                call write_xyz(v, angs_per_ps, comment, element, unit=300)                                    
            end if
        end do
        
    else if (integration == 'velocity-verlet') then
        allocate(fm(n_part,3))
        do while(t <= tf)
            call velocity_verlet(del_t)
            t = t + del_t
            nstep = nstep + 1
            if ( modulo(nstep,print_nstep) .eq. 0 ) then
                call calculate_potential(en)
                ke =  0.5d0*mass*sum(v**2)
                write(20,'(f16.8, 4f32.16)') t*time_ps, en*eps_kjpm, ke*eps_kjpm, (en+ke)*eps_kjpm
                ! writing initial coordinates to xyzfile.
                write(comment, '(a,f16.8)') "Time t (ps) = ", t*time_ps
                call write_xyz(x, sigma_angs, comment, element, unit=200)
                call write_xyz(v, angs_per_ps, comment, element, unit=300)                                    
            end if
        end do
    else 
        print*, "!!ERROR: Could not recognise the given integration integration", integration
    end if
    
    close(20)
    close(200)
    close(300)
        
    
    
end program main