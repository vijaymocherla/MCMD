program main
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    use io_module
    implicit none
    
    integer :: n_part
    integer :: nstep
    integer :: rc, t, tmax, tau, ci

    real(dp), allocatable :: v(:,:), vi(:,:), vj(:,:)
    real(dp), allocatable :: vt(:,:,:)
    real(dp), allocatable :: vacf(:)
    
    real(dp) ::  sigma_angs
    real(dp) ::  time_ps
    real(dp) ::  angs_per_ps 
    real(dp) ::  box(3)
    real(dp) ::  memfp
    real(dp) ::  del_t

    character (len=1024) :: comment
    character (len=1024) :: velocity_file
    character (len=1024) :: outfile
    character (len=1024) :: arg
    character (len=32)   :: boxvar
    character (len=2)    :: element

    sigma_angs = 3.405d0
    time_ps = 2.180d0
    angs_per_ps = sigma_angs/time_ps

    ci = 1
    do
        call get_command_argument(ci, arg) 
        if (trim(arg) == "-vel") then               ! .vel filename
            call get_command_argument(ci+1, arg) 
            velocity_file = trim(arg)
            ci = ci + 2
        else if (trim(arg) == "-o") then            ! .corr output filename
            call get_command_argument(ci+1, arg) 
            outfile = trim(arg)
            ci = ci + 2   
        else if (trim(arg) == "-nstep") then        ! input conversion factor    
            call get_command_argument(ci+1, arg) 
            read(arg,'(i16)') nstep
            ci = ci + 2 
        else if (trim(arg) == "-del_t") then        ! input conversion factor    
            call get_command_argument(ci+1, arg) 
            read(arg,'(f32.16)') del_t
            ci = ci + 2  
        else if (trim(arg) == "-sigma_angs") then   ! input conversion factor    
            call get_command_argument(ci+1, arg) 
            read(arg,'(f32.16)') sigma_angs
            ci = ci + 2   
        else
            exit
        end if
    end do
    open(100, file=velocity_file)
    call read_xyz(v, angs_per_ps, comment, element, 100, rc)
    read(comment,*) boxvar, box
    n_part = size(v, dim=1)
    ! calculate memory footprint
    memfp = real(nstep,8)*real(n_part,8)*3.0d0*8.0d0/1000000.0d0
    if (memfp > 2048.0d0)then
        print*, "Memory foot print exceeds set memory limi of 2Gb."
        stop
    end if
    ! intiating and allocaing arrays
    tmax = int(nstep/2)
    allocate(vt(n_part, 3, 0:nstep-1))
    allocate(vi(n_part, 3), vj(n_part, 3))
    allocate(vacf(0:tmax-1))
    vt = 0.0d0
    vacf = 0.0d0
    vt(:,:,0) = v(:,:)
    do t=1,nstep-1
        call read_xyz(v, angs_per_ps, comment, element, 100, rc)
        vt(:,:,t) = v(:,:)
    end do
    close(100)
    ! calculating autocorrelation function
    do tau=0,tmax-1
        do t=1,tmax-1
            vacf(tau) = vacf(tau) + sum(vt(:,:,t)*vt(:,:,t+tau))
        end do
    end do

    ! normalising vacf
    vacf = vacf/real(tmax,8)/real(n_part,8)

    ! write vacf to file
    open(200, file=outfile)
    del_t = time_ps*del_t
    write(200,'(i16)') tmax 
    write(200,'(2a32)') 't (ps)', 'C(t)'
    do t=0,tmax-1
        write(200,'(2f32.16)') del_t*t, vacf(t)
    end do
    close(200)
    
end program main