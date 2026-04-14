! lap2d_mod.f90
! 2D Laplace SLP/DLP evaluation for the exterior disk BIE.
!
! Representation formula (one target zt outside the disk):
!   u(zt) = (1/2pi) * sum_j [ ub(j)*Im(wzp(j)/(zt-z(j))) - unb(j)*(-log|z(j)-zt|)*awzp(j) ]
!
! Called for:
!   - List1: full np-node sum over all panel nodes
!   - List2: per-panel p-node sum for far panels
!
! Public routines:
!   lap2d_eval_r64 (np, z, wzp, awzp, ub, unb, zt, u)
!   lap2d_eval_r128(np, z, wzp, awzp, ub, unb, zt, u)

module lap2d_mod
  use utils_mod, only: r64, r128
  implicit none

contains

  subroutine lap2d_eval_r64(np, z, wzp, awzp, ub, unb, zt, u)
    integer(8),   intent(in)  :: np
    complex(r64), intent(in)  :: z(np), wzp(np), zt
    real(r64),    intent(in)  :: awzp(np), ub(np), unb(np)
    real(r64),    intent(out) :: u
    ! TODO: implement
  end subroutine lap2d_eval_r64

  subroutine lap2d_eval_r128(np, z, wzp, awzp, ub, unb, zt, u)
    integer(8),    intent(in)  :: np
    complex(r128), intent(in)  :: z(np), wzp(np), zt
    real(r128),    intent(in)  :: awzp(np), ub(np), unb(np)
    real(r128),    intent(out) :: u
    ! TODO: implement
  end subroutine lap2d_eval_r128

end module lap2d_mod
