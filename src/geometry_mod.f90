! geometry_mod.f90
! Panel discretization of 2D curves for the Laplace BIE.
! Currently supports: unit circle (z(s) = exp(2*pi*i*s), s in [0,1)).
! Outward normal for exterior domain: nz = -i*zp/|zp|.
!
! Public routines:
!   zinit_r64 (npan, p, sinter, sinterdiff, T, W, z, zp, zpp, nz, w, wzp, awzp, zinter)
!   zinit_r128(npan, p, sinter, sinterdiff, T, W, z, zp, zpp, nz, w, wzp, awzp, zinter)

module geometry_mod
  use utils_mod, only: r64, r128
  implicit none

contains

  subroutine zinit_r64(npan, p, sinter, sinterdiff, T, Wq, &
      z, zp, zpp, nz, w, wzp, awzp, zinter)
    integer(8),   intent(in)  :: npan, p
    real(r64),    intent(in)  :: sinter(npan+1), sinterdiff(npan), T(p), Wq(p)
    complex(r64), intent(out) :: z(npan*p), zp(npan*p), zpp(npan*p)
    complex(r64), intent(out) :: nz(npan*p), wzp(npan*p), zinter(npan+1)
    real(r64),    intent(out) :: w(npan*p), awzp(npan*p)
    real(r64)    :: s(npan*p), sdif, pi
    complex(r64) :: ci, ci2pi
    integer(8)   :: k, i1, i2

    pi    = acos(-1.0_r64)
    ci    = cmplx(0.0_r64, 1.0_r64, r64)
    ci2pi = 2.0_r64 * pi * ci

    ! Build global quadrature nodes s and weights w panel by panel
    do k = 1, npan
      i1 = (k-1)*p + 1
      i2 = k*p
      sdif = sinterdiff(k) / 2.0_r64
      s(i1:i2) = (sinter(k) + sinter(k+1)) / 2.0_r64 + sdif * T
      w(i1:i2) = sdif * Wq
    end do
    ! z(s) = exp(2*pi*i*s)
    z = exp(ci2pi * s)
    ! z'(s) = 2*pi*i*exp(2*pi*i*s)
    zp = ci2pi * z
    ! z''(s) = -(2*pi)^2 * exp(2*pi*i*s)
    zpp = -(2.0_r64*pi)**2 * z
    ! z at panel interfaces
    zinter = exp(ci2pi * sinter)
    ! outward normal for exterior domain
    nz = -ci * zp / abs(zp)
    ! weighted derivatives
    wzp  = w * zp
    awzp = abs(wzp)
  end subroutine zinit_r64

  subroutine zinit_r128(npan, p, sinter, sinterdiff, T, Wq, &
      z, zp, zpp, nz, w, wzp, awzp, zinter)
    integer(8),    intent(in)  :: npan, p
    real(r128),    intent(in)  :: sinter(npan+1), sinterdiff(npan), T(p), Wq(p)
    complex(r128), intent(out) :: z(npan*p), zp(npan*p), zpp(npan*p)
    complex(r128), intent(out) :: nz(npan*p), wzp(npan*p), zinter(npan+1)
    real(r128),    intent(out) :: w(npan*p), awzp(npan*p)
    real(r128)    :: s(npan*p), sdif, pi
    complex(r128) :: ci, ci2pi
    integer(8)    :: k, i1, i2

    pi    = acos(-1.0_r128)
    ci    = cmplx(0.0_r128, 1.0_r128, r128)
    ci2pi = 2.0_r128 * pi * ci

    do k = 1, npan
      i1 = (k-1)*p + 1
      i2 = k*p
      sdif = sinterdiff(k) / 2.0_r128
      s(i1:i2) = (sinter(k) + sinter(k+1)) / 2.0_r128 + sdif * T
      w(i1:i2) = sdif * Wq
    end do

    z = exp(ci2pi * s)
    zp = ci2pi * z
    zpp = -(2.0_r128*pi)**2 * z
    zinter = exp(ci2pi * sinter)
    nz = -ci * zp / abs(zp)
    wzp  = w * zp
    awzp = abs(wzp)
  end subroutine zinit_r128

end module geometry_mod
