module utils
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = real64

  private :: dp
  
contains
  
    subroutine linspace(start, end, arr)
    real(dp), intent(in) :: start
    real(dp), intent(in) :: end

    real(dp), dimension(:), intent(out) :: arr

    integer :: n 
    integer :: i

    n = size(arr)
    if (n==0) return
    if (n==1) then
       arr(1) = start
    end if

    do i=1,n
       arr(i) = start + (end - start) * (i-1)/(n-1)
    end do        
    
  end subroutine linspace

  subroutine logspace(start, end, arr)
    real(dp), intent(in) :: start, end
    real(dp), intent(inout)  :: arr(:)

    integer :: n
    integer :: i

    n = size(arr)
    call linspace(start, end, arr)

    do i=1,n
       arr(i) = 10**(arr(i))
    end do
    
  end subroutine logspace

  pure subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

    !yval = 10._dp**( ly1 + (lxval - lx1)*(ly2 - ly1)*norm )
  end subroutine linear_log_interp

  pure subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none
    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = ( y1*(x2 - xval) + y2*(xval - x1) ) * norm

  end subroutine linear_interp

  subroutine lin_interp(x, data, x_new, data_out)
    real(dp), intent(in) :: x(:), data(:)
    real(dp), intent(in) :: x_new
    real(dp), intent(out) :: data_out

    integer :: i, ii, itemp(1) 

    itemp = minloc(abs(x - x_new))
    ii = itemp(1)

    if (ii .eq. size(data)) then
       i = ii-1
       if (i.eq.0) print*, '1'
    else if (i .eq. 1 ) then
       i = ii
       if (i.eq.0) print*, '2'
    else if (x_new > x(ii)) then
       i = ii
       if (i.eq.0) print*, '3'
    else
       i = ii-1
       if (i.eq.0) print*, '4', ii, x_new
    endif

    if (i .eq. 0) then
       write(*,*) 'i IS ZERO'
    endif
    
   
    data_out = data(i) + (x_new - x(i))/(x(i+1) - x(i)) * (data(i+1) - data(i))
  end subroutine lin_interp
  
  subroutine bilinear_interp(x, y, data, x_new, y_new, data_out)
    real(dp), dimension(:), intent(in) :: x,y
    real(dp), intent(in)  :: x_new, y_new
    real(dp), dimension(:,:), intent(in) :: data
    
    real(dp), intent(out) :: data_out

    integer :: i,j, ii,jj
    integer :: itemp(1), jtemp(1)
    real(dp) :: v1(2), v2(2), v3(2), mat(2,2)
    
    ! Find nearest neighbour points

    !write(*,*) "SHAPE ", shape(data)
    
    itemp = minloc(abs(x - x_new))
    i = itemp(1)
    
    if (i .eq. size(data,1)) then
       ! Accounts for if x(i) > x(nx)
       ii = i-1
    else if (i .eq. 1) then
       ! Accounts for if x(i) < x(1)
       ii = i
    else if (x_new > x(i)) then
       ! Ensures x(i) < x_new < x(i+1)
       ii = i
    else
       ii = i-1
    endif
    
    jtemp = minloc(abs(y - y_new))
    j = jtemp(1)
    
    if (j .eq. size(data,2)) then
       jj = j-1
    else if (j .eq. 1) then
       jj = j
    else if (y_new > y(j)) then
       jj = j
    else
       jj = j-1
    endif

    v1 = (/ (x(ii+1)-x_new), (x_new-x(ii)) /)
    v2 = (/ (y(jj+1)-y_new), (y_new-y(jj)) /)
    mat = reshape( (/ data(ii,jj),  data(ii,jj+1),&
                     data(ii+1,jj), data(ii+1,jj+1) /), (/2,2/), &
                     order=(/2,1/) )

    v3 = matmul(mat, v2)
    data_out = dot_product(v1, v3)/(x(ii+1) - x(ii))/(y(jj+1) - y(jj))
       
  end subroutine bilinear_interp

!  subroutine linear_log_extrap(xval, x1, x2, y1, y2, yval)

!    real(dp), intent(in) :: xval, y1, y2, x


  subroutine writeline
    write(*,*) '-----------------------------------------------------'
  end subroutine writeline
  
end module utils
