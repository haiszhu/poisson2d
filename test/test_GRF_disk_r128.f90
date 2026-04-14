! test_GRF_disk_r128.f90
! Pure Fortran quad-precision GRF disk test.
! Mirrors test_GRF_disk_r64.f90 with r128 types/constants.
! Residual metric: max|u| on outside targets, where u should be near zero.

program test_GRF_disk_r128
  use utils_mod
  use geometry_mod
  use specialquad_mod
  implicit none

  integer(8), parameter :: npan=60, p=32, np=npan*p
  integer(8), parameter :: Ng=100, Ngtot=Ng*Ng
  real(r128), parameter :: dlim=0.8_r128, coef1=0.70_r128, coef2=1.00_r128
  integer(8), parameter :: pow1=8, pow2=7

  ! Panel arrays
  real(r128)    :: sinter(npan+1), sinterdiff(npan), T(p), W(p), Dgl(p,p)
  real(r128)    :: wq(np), awzp(np), ub(np), unb(np), panlen(npan)
  complex(r128) :: z(np), zp(np), zpp(np), nz(np), wzp(np), zinter(npan+1)

  ! Target arrays
  complex(r128) :: zg(Ngtot)
  real(r128)    :: xg(Ng), yg(Ng), u(Ngtot)
  integer(8)    :: mylist(Ngtot), list1(Ngtot), list2(Ngtot), nl1, nl2, nout

  ! Work
  integer(8)    :: k, kk, j, jj, i1, i2
  real(r128)    :: d, err1, err2, unum, unumS, unumK, pi, fx, fy
  complex(r128) :: zt(1), tmp1, tmp2, ci
  real(r128)    :: As(p,1), A1(p,1), A2(p,1), A3(p,1), A4(p,1)
  complex(r128) :: Ad(p,1)
  logical       :: is_out

  pi = acos(-1.0_r128)
  ci = cmplx(0.0_r128, 1.0_r128, r128)

  ! ---- panel setup ----
  call gauss_r128(p, T, W, Dgl)
  do k = 1, npan+1
    sinter(k) = real(k-1, r128) / real(npan, r128)
  end do
  do k = 1, npan
    sinterdiff(k) = 1.0_r128 / real(npan, r128)
  end do
  call zinit_r128(npan, p, sinter, sinterdiff, T, W, z, zp, zpp, nz, wq, wzp, awzp, zinter)

  ! ---- boundary data ----
  do j = 1, np
    ub(j) = coef1*real(z(j)**pow1, r128) + coef2*aimag(z(j)**pow2)
    tmp1  = real(pow1, r128) * z(j)**(pow1-1)
    tmp2  = real(pow2, r128) * z(j)**(pow2-1)
    fx    = coef1*real(tmp1, r128) + coef2*aimag(tmp2)
    fy    = coef1*real(ci*tmp1, r128) + coef2*aimag(ci*tmp2)
    unb(j) = real(nz(j), r128)*fx + aimag(nz(j))*fy
  end do

  ! ---- panel arc lengths ----
  do kk = 1, npan
    i1 = (kk-1)*p + 1
    i2 = kk*p
    panlen(kk) = sum(awzp(i1:i2))
  end do

  ! ---- target grid ----
  do k = 1, Ng
    xg(k) = -1.2_r128 + 2.4_r128*real(k-1, r128)/real(Ng-1, r128)
    yg(k) = xg(k)
  end do
  do k = 1, Ng
    do j = 1, Ng
      zg((k-1)*Ng + j) = cmplx(xg(k), yg(j), r128)
    end do
  end do

  ! ---- classify targets ----
  nl1 = 0
  nl2 = 0
  nout = 0
  do k = 1, Ngtot
    mylist(k) = 0
    is_out = abs(zg(k)) > 1.0_r128
    if (.not. is_out) cycle
    nout = nout + 1
    mylist(k) = 1
    do kk = 1, npan
      i1 = (kk-1)*p + 1
      i2 = kk*p
      d = minval(abs(z(i1:i2) - zg(k))) / panlen(kk)
      if (d < dlim) then
        mylist(k) = 2
        exit
      end if
    end do
    if (mylist(k) == 1) then
      nl1 = nl1 + 1
      list1(nl1) = k
    else
      nl2 = nl2 + 1
      list2(nl2) = k
    end if
  end do

  ! ---- list1: full direct sum ----
  do jj = 1, nl1
    k = list1(jj)
    unum = 0.0_r128
    do j = 1, np
      unumS = -log(abs(z(j)-zg(k))) * awzp(j) / (2.0_r128*pi)
      unumK =  aimag(wzp(j)/(zg(k)-z(j))) / (2.0_r128*pi)
      unum = unum + ub(j)*unumK - unb(j)*unumS
    end do
    u(k) = abs(unum)
  end do

  ! ---- list2: far direct + near special quad ----
  do jj = 1, nl2
    k = list2(jj)
    unum = 0.0_r128
    do kk = 1, npan
      i1 = (kk-1)*p + 1
      i2 = kk*p
      d = minval(abs(z(i1:i2)-zg(k))) / panlen(kk)

      if (d > dlim) then
        do j = i1, i2
          unumS = -log(abs(z(j)-zg(k))) * awzp(j) / (2.0_r128*pi)
          unumK =  aimag(wzp(j)/(zg(k)-z(j))) / (2.0_r128*pi)
          unum = unum + ub(j)*unumK - unb(j)*unumS
        end do
      else
        zt(1) = zg(k)
        call sdspecialquad_r128(1_8, zt, p, z(i1:i2), nz(i1:i2), wzp(i1:i2), &
             zinter(kk), zinter(kk+1), 1_8, As, Ad, A1, A2, A3, A4)
        do j = 1, p
          unum = unum + ub(i1+j-1)*real(Ad(j,1), r128) - unb(i1+j-1)*As(j,1)
        end do
      end if
    end do
    u(k) = abs(unum)
  end do

  ! ---- GRF residual diagnostics (outside only): u should be near zero ----
  err1 = 0.0_r128
  err2 = 0.0_r128
  do jj = 1, nl1
    k = list1(jj)
    err1 = max(err1, u(k))
  end do
  do jj = 1, nl2
    k = list2(jj)
    err2 = max(err2, u(k))
  end do

  write(*,'(A,I0,A,I0,A,I0)') 'Counts: outside=', nout, ', list1=', nl1, ', list2=', nl2
  write(*,'(A,1PE12.4,A,0PF8.2)') 'List1 GRF residual max|u|: ', err1, '  log10=', log10(err1 + epsilon(1.0_r128))
  write(*,'(A,1PE12.4,A,0PF8.2)') 'List2 GRF residual max|u|: ', err2, '  log10=', log10(err2 + epsilon(1.0_r128))

  if (nl1 == 0 .or. nl2 == 0) then
    write(*,*) 'test_GRF_disk_r128: invalid target classification'
    stop 2
  end if

  write(*,*) 'test_GRF_disk_r128: completed'

end program test_GRF_disk_r128
