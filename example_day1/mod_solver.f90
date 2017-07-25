module mod_solver
  use mod_types, only: SI, DP
  use mod_rhs, only: func
  implicit none

contains

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

end module mod_solver
