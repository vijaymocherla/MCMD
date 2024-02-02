module io_module
    ! To Do:
    !   - integrate all message passing to out files
    !   - simplify interface for .out files
    use, intrinsic :: iso_fortran_env, only:dp=>real64
    implicit none
    private

    public :: read_xyz, write_xyz

    !public data
    
    contains 

    subroutine read_xyz(array, alpha, comment, element, unit, rc)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        integer, intent(in)  :: unit    ! i/o unit for file
        real(dp), intent(in) :: alpha   ! conversion factor b/w reduced units and normal units
        integer, intent(out) :: rc      ! iostat variable
        character (len=2), intent(out) :: element
        character (len=1024), intent(out) :: comment 
        real(dp), allocatable, intent(out) :: array(:,:)
        integer :: n, i
        read(unit, *, iostat=rc) n
        if (.not. is_iostat_end(rc)) then
            allocate(array(n,3))
            read(unit, '(a)') comment
            do i=1,n
                read(unit, fmt='(a2,3f32.16)', iostat=rc) element, array(i, :)
            end do
            array = array/alpha
        end if 
    end subroutine

    subroutine write_xyz(array, alpha, comment, element, unit)
        use, intrinsic :: iso_fortran_env, only:dp=>real64
        implicit none
        integer, intent(in)  :: unit    ! i/o unit for file
        real(dp), intent(in) :: alpha   ! conversion factor b/w reduced units and normal units
        character (len=2), intent(in) :: element
        character (len=1024), intent(in) :: comment 
        real(dp), allocatable, intent(in) :: array(:,:)
        integer :: n, i

        n = size(array, dim=1)
        write(unit, *) n
        write(unit, '(a)') trim(comment)
        do i=1,n
            write(unit,'(a2,3f32.16)') element, alpha*array(i, :)
        end do
    end subroutine

end module io_module