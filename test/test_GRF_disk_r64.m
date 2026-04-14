% GRF_disk_r64.m
% Double-precision GRF disk test via poisson2d MEX interface.
% Mirrors GRF_disk.m; uses Fortran for gauss, zinit, lap2d_eval.
% Near-panel special quadrature falls back to MATLAB SDspecialquad (from GRF_disk.m).
%
% Requires: poisson2d_mex built (make mex), SDspecialquad.m in this folder.

testdir = fileparts(mfilename('fullpath'));
addpath(fullfile(testdir, '../matlab'))

p     = 16;
npan  = 60;
dlim  = 0.8;
coef1 = 0.70;  pow1 = 8;
coef2 = 1.00;  pow2 = 7;

% --- panel setup via MEX ---
[T, W, ~] = gauss_mex(p);
sinter     = linspace(0,1,npan+1)';
sinterdiff = ones(npan,1)/npan;
[z, zp, zpp, nz, w, wzp, awzp, zinter] = zinit_mex(npan, p, sinter, sinterdiff, T, W);
% [z2,zp2,~,nz2,~,wzp2,awzp2,zinter2] = zinit(sinter,sinterdiff,T,W,npan);

% target grid
Ngtot = Ng*Ng;
xg = linspace(-1.2,1.2,Ng);
yg = linspace(-1.2,1.2,Ng);
zg = zeros(Ngtot,1);
for k = 1:Ng
  zg((k-1)*Ng+(1:Ng)) = xg(k) + 1i*yg;
end

% boundary data from high-order harmonic test function (ReviewerX_Q4Q5 style)
f  = coef1*real(z.^pow1) + coef2*imag(z.^pow2);
fx = coef1*real(pow1*z.^(pow1-1)) + coef2*imag(pow2*z.^(pow2-1));
fy = coef1*real(1i*pow1*z.^(pow1-1)) + coef2*imag(1i*pow2*z.^(pow2-1));
ub  = f;
unb = real(nz).*fx + imag(nz).*fy;

% classify outside targets into list1/list2
panlen = zeros(npan,1);
for kk = 1:npan
  myind = (kk-1)*p + (1:p);
  panlen(kk) = sum(awzp(myind));
end

inout = Ineval(zg);
mylist = inout;
for k = 1:Ngtot
  if inout(k)
    for kk = 1:npan
      myind = (kk-1)*p + (1:p);
      d = min(abs(z(myind)-zg(k)))/panlen(kk);
      if d < dlim
        mylist(k) = 2;
      end
    end
  end
end
mylist1 = find(mylist==1)';
mylist2 = find(mylist==2)';

%
u = eps*ones(Ngtot,1);

%
disp('List1 starts')
for k = mylist1
  UnumS = -log(abs(z-zg(k))).*awzp/2/pi;
  UnumK = imag(wzp./(zg(k)-z))/2/pi;
  Unum = ub.'*UnumK - unb.'*UnumS;
  u(k) = abs(Unum);
end

disp('List2 starts')
for k = mylist2
  Unum = 0;
  for kk = 1:npan
    myind = (kk-1)*p + (1:p);
    zsc_kk = z(myind);
    wzp_kk = wzp(myind);
    awzp_kk = abs(wzp_kk);
    ub_kk = ub(myind);
    unb_kk = unb(myind);
    d = min(abs(zsc_kk-zg(k)))/panlen(kk);

    if d > dlim
      UnumS = -log(abs(zsc_kk-zg(k))).*awzp_kk/2/pi;
      UnumK = imag(wzp_kk./(zg(k)-zsc_kk))/2/pi;
      Unum = Unum + ub_kk.'*UnumK - unb_kk.'*UnumS;
    else
      a_kk = zinter(kk);
      b_kk = zinter(kk+1);
      nz_kk = nz(myind);

      % tt = struct('x', zg(k));
      % sk = struct('x', zsc_kk, 'nx', nz_kk, 'wxp', wzp_kk);
      % [Sclose, Dclose] = SDspecialquad(tt, sk, a_kk, b_kk, 'e');
      % Unum = Unum + ub_kk.'*real(Dclose)' - unb_kk.'*real(Sclose)';

      [Sclose, Dclose] = sdspecialquad_mex( 1, zg(k), p, zsc_kk, nz_kk, wzp_kk, a_kk, b_kk, 1);
      Unum = Unum + ub_kk.'*real(Dclose(:,1)) - unb_kk.'*Sclose(:,1);
    end
  end
  u(k) = abs(Unum);
end

u(u<eps)=eps;
F1=zeros(Ng);
for k=1:Ng
  F1(1:Ng,k)=log10(u((k-1)*Ng+(1:Ng)));
end
figure(1),clf,
imagesc(xg,yg,F1);      
colormap(flipud(pink))
axis xy
colorbar
hold on

keyboard

function iout = Ineval(zg)
iout = zeros(numel(zg),1);
iout(abs(zg) > 1) = 1;
end
