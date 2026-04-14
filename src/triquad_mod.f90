! triquad_mod.f90
! Vioreanu-Rokhlin quadrature nodes and weights on the reference triangle
! (0,0)-(1,0)-(0,1), compiled in from koorn-uvs-dat.txt / koorn-wts-dat.txt.
!
! Public:
!   koorn_uvs_wts_r64(norder, npols, uvs, wts)
!     norder       : polynomial order (0..20)
!     npols        : (norder+1)*(norder+2)/2  [output]
!     uvs(2,npols) : (u,v) coordinates on reference triangle  [output]
!     wts(npols)   : integration weights (sum = 0.5)          [output]

module triquad_mod
  use utils_mod, only: r64
  implicit none
  private
  public :: koorn_uvs_wts_r64

contains

  subroutine koorn_uvs_wts_r64(norder, npols, uvs, wts)
    integer(8), intent(in)  :: norder
    integer(8), intent(out) :: npols
    real(r64),  intent(out) :: uvs(2,*), wts(*)
    npols = (norder+1_8)*(norder+2_8)/2_8
    include 'koorn-uvs-dat.txt'
    include 'koorn-wts-dat.txt'
  end subroutine koorn_uvs_wts_r64

end module triquad_mod
