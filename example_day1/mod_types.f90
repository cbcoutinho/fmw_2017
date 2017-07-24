module mod_types
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none

  private
  public :: SI, DI, SP, DP

  integer, parameter :: SI = int32
  integer, parameter :: DI = int64
  integer, parameter :: SP = real32
  integer, parameter :: DP = real64

contains

end module mod_types
