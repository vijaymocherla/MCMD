program main
    
    use, intrinsic :: iso_fortran_env, only:dp=>real64    
    use io_module
    use montecarlo
    
    implicit none
        
    real(dp) :: en
    logical  :: decision
    
    integer :: i, icycle, rc
    real(dp) :: rc2i, rc6i
    real(dp) :: sigma_angs, eps_kjpm ! conversion factors

    
    ! Input parameters types
    integer  :: ncycle                        ! no. of Monte Carlo steps    
    integer  :: print_nstep                   ! print every n-steps
    character(len=1024) :: initial_positions  ! initial configuration of atoms
    character(len=1024) :: outfile            ! output file to print energies  
    character(len=1024) :: xyz_file           ! xyz file to store trajectory 
    character(len=1024) :: comment            ! comment line for xyz file 
    character(len=2)    :: element            ! atom type 

    integer  :: nseed
    integer  :: moves
    character (len=32) :: boxvar 

    ! input from standard namelist
    namelist /mc_params/ initial_positions, outfile, xyz_file, &
                         temperature, delta_max, ncycle, print_nstep, &
                         sigma_angs, eps_kjpm, r_cut
                         
    open(10, file='mc_lj.inp')
        read(10, nml=mc_params, iostat=rc)
        if (rc /= 0) then
            print*, '!! I/O ERROR: Can not read input file. IOSTAT=', rc
            stop
        end if
    close(10)
    
    ! Setting up the random seed
    nseed = 1000
    call random_seed(size=nseed) 

    ! read initial configuration
    open(100, file=initial_positions)
        call read_xyz(x, sigma_angs, comment, element, unit=100, rc=rc)
        read(comment,'(a32,3f16.8)') boxvar, box
        box = box/sigma_angs   
    close(100)
    n_part = size(x, dim=1)
    ! initialising variables
    en = 0.0d0
    ! calculate cutoff energy
    rc2i = 1.0d0 / r_cut**2
    rc6i = rc2i**3
    e_cut = 4.0d0*rc6i*(rc6i - 1.0d0)
    moves = 0.0d0
    
    ! start mc
    open(20,file=outfile)
    open(200,file=xyz_file)
    call calculate_energy(en)
    write(20,*) n_part, ncycle, en*eps_kjpm
    write(20,*) "! comment line"
    write(20,'(a8, a8, a12, a32)') "icycle", "i", "decision", "Potential Energy (in KJ/mol)" 
    write(comment, '(a32,3f16.8)') "'!  box_size (in Angstrom)  =  '", box*sigma_angs
    call write_xyz(x, sigma_angs, comment, 'Ar', unit=200)
    do icycle = 1, ncycle
        call mcmove(i, decision, en)
        write(20,'(i8, i8, l12, f32.16)') icycle, i, decision, en*eps_kjpm
        if (decision) then
            moves = moves + 1
        end if
        if (modulo(icycle,print_nstep) .eq. 0) then
            write(comment,*) 'icycle = ', icycle
            call write_xyz(x, sigma_angs, comment, 'Ar', unit=200)
        end if
    end do
    close(200)
    write(20,*)
    write(20,*) "Final energy (in KJ/mol): ", en
    write(20,*) "Accepted moves: ", moves
    write(20,*) "Ratio: ", real(moves)/real(ncycle)
    close(20)

end program main
