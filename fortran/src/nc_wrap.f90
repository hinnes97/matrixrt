module nc_wrap_mod
  use netcdf
  use io,     only : handle_err
  use params, only : dp
  implicit none

  interface get_var
     module procedure get_var_1d
     module procedure get_var_2d
     module procedure get_var_3d
     module procedure get_var_4d
  end interface
  
contains
  
  subroutine open_file(file_name, ncid)
    character(*), intent(in) :: file_name
    integer,     intent(out) :: ncid

    integer :: status

    status = nf90_open(file_name, 0, ncid)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine open_file

  subroutine close_file(ncid)
    integer, intent(in) :: ncid

    integer :: status

    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine close_file

  subroutine get_dimension_length(ncid, dim, N)
    integer, intent(in)      :: ncid
    character(*), intent(in) :: dim
    integer, intent(out)     :: N

    integer :: status, dimid
    
    status = nf90_inq_dimid(ncid, dim, dimid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimid, len=N)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine get_dimension_length

  subroutine get_var_1d(ncid, var, data)
    integer, intent(in) :: ncid
    character(*), intent(in) :: var
    real(dp), allocatable, intent(out) :: data(:)

    integer :: ndim, status, varid
    integer :: lens(nf90_max_var_dims)
    
    call get_var_info(ncid, var, varid, ndim, lens)

    if (ndim .ne. 1) then
       stop 'Variable has wrong rank - must be 1D'
    else
       allocate( data(lens(1)) )
    endif
        
    status = nf90_get_var(ncid, varid, data)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine get_var_1d

  subroutine get_var_2d(ncid, var,data)
    integer, intent(in) :: ncid
    character(*), intent(in) :: var
    real(dp), allocatable, intent(out) :: data(:,:)

    integer :: ndim, status, varid
    integer :: lens(nf90_max_var_dims)
    
    call get_var_info(ncid, var, varid, ndim, lens)

    if (ndim .ne. 2) then
       stop 'Variable has wrong rank - must be 2D'
    else
       allocate( data(lens(1), lens(2)) )
    endif

    status = nf90_get_var(ncid, varid, data)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine get_var_2d

  subroutine get_var_3d(ncid, var,data)
    integer, intent(in) :: ncid
    character(*), intent(in) :: var
    real(dp), allocatable, intent(out) :: data(:,:,:)

    integer :: ndim, status, varid
    integer :: lens(nf90_max_var_dims)
    
    call get_var_info(ncid, var, varid, ndim, lens)

    if (ndim .ne. 3) then
       stop 'Variable has wrong rank - must be 3D'
    else
       allocate( data(lens(1), lens(2), lens(3)) )
    endif
    
    status = nf90_get_var(ncid, varid, data)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine get_var_3d

  subroutine get_var_4d(ncid, var, data)
    integer, intent(in) :: ncid
    character(*), intent(in) :: var
    real(dp), allocatable, intent(out) :: data(:,:,:,:)

    integer :: ndim, status, varid
    integer :: lens(nf90_max_var_dims)
    
    call get_var_info(ncid, var, varid, ndim, lens)

    if (ndim .ne. 4) then
       stop 'Variable has wrong rank - must be 4D'
    else
       allocate( data(lens(1), lens(2), lens(3), lens(4)) )
    endif
        
    status = nf90_get_var(ncid, varid, data)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine get_var_4d

  subroutine get_var_info(ncid, var, varid, ndim, lens)
    integer, intent(in) :: ncid
    character(*), intent(in) :: var
    integer, intent(out) :: varid, ndim
    integer, dimension(nf90_max_var_dims) :: lens

    integer, dimension(nf90_max_var_dims) :: dimids
    integer :: n, status
    
    status = nf90_inq_varid(ncid, var, varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_variable(ncid, varid, ndims=ndim, dimids=dimids)
    if (status /= nf90_noerr) call handle_err(status)

    do n=1,ndim
       status = nf90_inquire_dimension(ncid, dimids(n), len=lens(n))
       if (status /= nf90_noerr) call handle_err(status)
    enddo

  end subroutine get_var_info
      
end module nc_wrap_mod
