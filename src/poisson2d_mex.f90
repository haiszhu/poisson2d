! poisson2d_mex.f90
! Standalone (non-module) MEX-facing wrappers.
!
! gfortran module procedures export as __mod_MOD_foo (module-mangled);
! mwrap expects plain symbol _foo_ (trailing underscore on macOS).
! Each wrapper imports the module procedure under a local USE alias.

! ------------------------------------------------------------------
! gauss
! ------------------------------------------------------------------
subroutine gauss_r64(n, tgl, wgl, Dgl)
  use utils_mod, only: g => gauss_r64
  implicit none
  integer(8), intent(in)    :: n
  real(8),    intent(inout) :: tgl(n), wgl(n), Dgl(n,n)
  call g(n, tgl, wgl, Dgl)
end subroutine gauss_r64

! ------------------------------------------------------------------
! zinit: panel discretization of the unit circle
! ------------------------------------------------------------------
subroutine zinit_r64(npan, p, sinter, sinterdiff, T, Wq, &
    z, zp, zpp, nz, w, wzp, awzp, zinter)
  use geometry_mod, only: zi => zinit_r64
  implicit none
  integer(8), intent(in)    :: npan, p
  real(8),    intent(in)    :: sinter(npan+1), sinterdiff(npan), T(p), Wq(p)
  complex(8), intent(inout) :: z(npan*p), zp(npan*p), zpp(npan*p)
  complex(8), intent(inout) :: nz(npan*p), wzp(npan*p), zinter(npan+1)
  real(8),    intent(inout) :: w(npan*p), awzp(npan*p)
  call zi(npan, p, sinter, sinterdiff, T, Wq, z, zp, zpp, nz, w, wzp, awzp, zinter)
end subroutine zinit_r64

! ------------------------------------------------------------------
! lap2d_eval: SLP/DLP sum at one target.
! Target zt split into (ztre, ztim) — mwrap does not support complex scalars.
! ------------------------------------------------------------------
subroutine lap2d_eval_r64(np, z, wzp, awzp, ub, unb, ztre, ztim, u)
  use lap2d_mod, only: ev => lap2d_eval_r64
  implicit none
  integer(8), intent(in)  :: np
  complex(8), intent(in)  :: z(np), wzp(np)
  real(8),    intent(in)  :: awzp(np), ub(np), unb(np), ztre, ztim
  real(8),    intent(out) :: u
  complex(8) :: zt
  zt = cmplx(ztre, ztim, 8)
  call ev(np, z, wzp, awzp, ub, unb, zt, u)
end subroutine lap2d_eval_r64

! ------------------------------------------------------------------
! sdspecialquad_r64: Helsing close-evaluation weights for nt targets, 1 panel.
!
! Outputs (p×nt, Fortran column-major):
!   As(p,nt)   : SLP weights (real)
!   A(p,nt)    : DLP weights (complex; take real for Laplace)
!   A1(p,nt)   : Re(Az), x-deriv
!   A2(p,nt)   : -Im(Az), y-deriv
!   A3(p,nt)   : Re(Azz)
!   A4(p,nt)   : -Im(Azz)
!
! Panel endpoints a, b: passed as dcomplex[1] from MATLAB.
! iside: 1 = exterior, 0 = interior.
! ------------------------------------------------------------------
subroutine sdspecialquad_r64(nt, zt, p, zsrc, nzsrc, wzp, a, b, iside, &
    As, Ad, A1, A2, A3, A4)
  use specialquad_mod, only: sq => sdspecialquad_r64
  implicit none
  integer(8), intent(in)    :: nt, p
  complex(8), intent(in)    :: zt(nt), zsrc(p), nzsrc(p), wzp(p)
  complex(8), intent(in)    :: a(1), b(1)
  integer(8), intent(in)    :: iside
  real(8),    intent(inout) :: As(p,nt), A1(p,nt), A2(p,nt), A3(p,nt), A4(p,nt)
  complex(8), intent(inout) :: Ad(p,nt)
  call sq(nt, zt, p, zsrc, nzsrc, wzp, a(1), b(1), iside, As, Ad, A1, A2, A3, A4)
end subroutine sdspecialquad_r64

! ------------------------------------------------------------------
! koorn_uvs_wts_r64: Vioreanu-Rokhlin nodes/weights on reference triangle.
! npols not passed -- MATLAB pre-allocates based on (norder+1)*(norder+2)/2.
! ------------------------------------------------------------------
subroutine koorn_uvs_wts_r64(norder, uvs, wts)
  use triquad_mod, only: kuw => koorn_uvs_wts_r64
  implicit none
  integer(8), intent(in)    :: norder
  real(8),    intent(inout) :: uvs(2,*), wts(*)
  integer(8) :: npols
  call kuw(norder, npols, uvs, wts)
end subroutine koorn_uvs_wts_r64
