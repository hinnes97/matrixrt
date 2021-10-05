! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
program cdl2cdf

  use realtype_rd
  use def_std_io_icf
  use filenamelength_mod, ONLY: filenamelength

  implicit none

  integer :: ierr = 0                ! Error flag
  integer :: i, j                    ! Loop variables
  character (LEN=filenamelength) :: cdl_name, cdf_name

  integer :: nd_cdl_dimen  = 7       ! Array size for netCDF dimensions
  integer :: nd_cdl_dimen_size = 768 ! Array size for number of values
                                     !  in a netCDF-dimension
  integer :: nd_cdl_data = 768*768   ! Array size for netCDF data
  integer :: nd_cdl_var =7           ! Number of netCDF variables

  integer :: n_dimension       ! Number of netCDF dimesions
  integer :: n_var             ! Number of netCDF variables

  character :: &
       dimension_name(7)*24 ,& ! Names of dimensions
       dimension_type(7)*6  ,& ! Types of dimensions
       dimension_long(7)*40 ,& ! Long names of dimensions
       dimension_unit(7)*20 ,& ! Units of dimensions
       var_name(7)*24       ,& ! Name of variable
       var_type(7)*6        ,& ! Type of variable
       var_long(7)*40       ,& ! Long name of variable
       var_unit(7)*20          ! Unit of variable

  ! Values of integral dimensions
  integer :: dimension_array_int(768, 7)

  ! Values of floating point dimensions
  real(RealK) :: dimension_array_fl(768, 7)

  ! Sizes of dimensions
  integer :: dimension_size(7)

  ! Number of data elements
  integer :: n_data(7)

  ! Number of dimensions used by variable
  integer :: n_dimension_var(7)

  ! List of dimensions used in the variable
  integer :: list_dimension_var(7, 7)

  ! Integral data fields
  integer :: data_int(768*768, 7)

  ! Floating point data fields
  real(RealK) :: data_fl(768*768, 7)


  ! Initialize the character strings
  do i=1, nd_cdl_dimen
     do j=1, len(dimension_type(i))
        dimension_type(i)(j:j)=' '
     enddo
  enddo
  do i=1, nd_cdl_var
     do j=1, len(var_type(i))
        var_type(i)(j:j)=' '
     enddo
  enddo

  write(iu_stdout, '(a)') 'Enter the input CDL filename:'
  read(iu_stdin, '(a)') cdl_name

  call read_cdl(ierr, cdl_name,                                       &
       nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,      &
       n_dimension, dimension_name, dimension_type, dimension_unit,   &
       dimension_long, dimension_size,                                &
       dimension_array_int, dimension_array_fl,                       &
       n_var, var_name, var_type, var_unit, var_long,                 &
       n_dimension_var, list_dimension_var,                           &
       n_data, data_int, data_fl )

  do i=1,n_dimension
     j=1
     do while ( j < len_trim(dimension_name(i)))
        if (dimension_name(i)(j:j) == ' ') dimension_name(i)(j:j) = '_'
        j=j+1
     end do
  end do

  do i=1,n_var
     j=1
     do while ( j < len_trim(var_name(i)))
        if (var_name(i)(j:j) == ' ') var_name(i)(j:j) = '_'
        j=j+1
     end do
  end do

  ! The list of dimensions will be in C order and need to be reversed
  ! before passing to the writing routine:
  do j = 1, n_var
    list_dimension_var(1:n_dimension_var(j), j) =                     &
      list_dimension_var(n_dimension_var(j):1:-1, j)
  end do

  write(iu_stdout, '(a)') 'Enter the output netCDF filename:'
  read(iu_stdin, '(a)') cdf_name

  call write_cdf(ierr, cdf_name,                                      &
       nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,      &
       n_dimension, dimension_name, dimension_type, dimension_unit,   &
       dimension_long, dimension_size,                                &
       dimension_array_int, dimension_array_fl,                       &
       n_var, var_name, var_type, var_unit, var_long,                 &
       n_dimension_var, list_dimension_var,                           &
       n_data, data_int, data_fl )

end program cdl2cdf
