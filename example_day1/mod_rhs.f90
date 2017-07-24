module mod_rhs
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

end module mod_rhs
