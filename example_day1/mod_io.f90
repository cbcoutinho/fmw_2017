module mod_io
  use mod_types
  use mod_rhs
  implicit none

contains

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

end module mod_io
