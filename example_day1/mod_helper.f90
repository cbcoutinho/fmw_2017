module mod_helper
  use mod_types
  implicit none

contains

  function func( j, x ) result ( d )

    integer(kind=SI), intent(in) :: j
    real(kind=DP), intent(in) :: x(:)
    real(kind=DP), allocatable :: d(:)

    allocate(d, source=x)
    d(:) = 0.0_DP

    return
  end function func

  subroutine fd1d_heat_explicit( x, t, dt, cfl, h, h_new )

    real(kind=DP), intent(in) :: cfl
    real(kind=DP), intent(in) :: dt
    real(kind=DP), intent(in) :: t
    real(kind=DP), intent(in) :: h(:)
    real(kind=DP), intent(in) :: x(:)
    real(kind=DP), intent(out) :: h_new(:)

    integer(kind=SI) :: j
    real(kind=DP), allocatable :: f(:)

    f = func( j, x )

    h_new(1) = 0.0_DP

    do j = 2, ubound(x,1) - 1
      h_new(j) = h(j) + dt * f(j) + cfl * ( h(j-1) - 2.0_DP * h(j) + h(j+1) )
    end do

    ! set the boundary conditions again
    h_new(1) = 90.0_DP
    h_new(ubound(h_new,1)) = 70.0_DP
  end subroutine fd1d_heat_explicit

  subroutine fd1d_heat_explicit_cfl( k, t_num, t_min, t_max, x_num, x_min, x_max, cfl )

    real(kind=DP), intent(in) :: k
    real(kind=DP), intent(in) :: t_max
    real(kind=DP), intent(in) :: t_min
    integer(kind=SI), intent(in) :: t_num
    real(kind=DP), intent(in) :: x_max
    real(kind=DP), intent(in) :: x_min
    integer(kind=SI), intent(in) :: x_num
    real(kind=DP), intent(out) :: cfl

    real(kind=DP) :: dt
    real(kind=DP) :: dx

    dx = ( x_max - x_min ) / dble( x_num - 1 )
    dt = ( t_max - t_min ) / dble( t_num - 1 )

    cfl = k * dt / dx / dx

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl

    return
  end subroutine fd1d_heat_explicit_cfl

  subroutine r8mat_write( output_filename, table )

    integer(kind=SI) :: j
    character (len=*) :: output_filename
    integer(kind=SI) :: output_unit_id
    character (len=30) :: string
    real(kind=DP) :: table(:,:)

    integer(kind=SI) :: m, n, err
    m = size(table,1)
    n = size(table,2)

    output_unit_id = 100
    open( unit = output_unit_id, file = output_filename, status = 'replace' , iostat=err)
    if ( err /= 0 ) then
      print *, 'Error opening output file, ' // trim(output_filename)
      print *, 'Error code: ', err
      stop
    endif

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)', iostat=err ) '(', m, 'g', 24, '.', 16, ')'
    if ( err /= 0 ) then
      print*, 'Error writing to string'
    endif

    do j = 1, n
      write ( output_unit_id, string, iostat=err) table(1:m, j)
      ! print*, table(1:m,j)
      if ( err /= 0 ) then
        print*, 'Error writing to string'
      endif
    end do

    close( unit = output_unit_id, iostat=err)
    if ( err /= 0 ) then
      print*, 'Close encountered an error closing output'
    end if
  end subroutine r8mat_write

  subroutine r8vec_linspace ( a_first, a_last, a )

    real(kind=DP) :: a(:)
    real(kind=DP) :: a_first
    real(kind=DP) :: a_last
    integer(kind=SI) :: i

    integer(kind=SI) :: n
    n = size(a, 1)

    do i = 1, n
      a(i) = ( dble( n - i ) * a_first + dble( i - 1 ) * a_last ) / dble( n - 1 )
    end do

    return
  end subroutine r8vec_linspace

  subroutine r8vec_write ( output_filename, x )

    integer(kind=SI) :: j
    character(len=*) :: output_filename
    integer(kind=SI) :: output_unit_id
    real(kind=DP) :: x(:)

    integer(kind=SI) :: n
    n = size(x,1)

    output_unit_id = 11
    open( unit = output_unit_id, file = output_filename, status = 'replace' )

    do j = 1, n
      write ( output_unit_id, '(2x,g24.16)' ) x(j)
    end do

    close ( unit = output_unit_id )

    return
  end subroutine r8vec_write

end module mod_helper
