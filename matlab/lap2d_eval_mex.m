function [u] = lap2d_eval_mex(np, z, wzp, awzp, ub, unb, zt, u)
np   = double(np);
ztre = double(real(zt));
ztim = double(imag(zt));
u    = double(0);
mex_id_ = 'lap2d_eval_r64(c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[x])';
[u] = poisson2d_mex(mex_id_, np, z, wzp, awzp, ub, unb, ztre, ztim, u, 1, np, np, np, np, np, 1, 1, 1);
end

% --------------------------------------------------------------------------
% sdspecialquad_mex: Helsing close-evaluation weights for nt targets, 1 panel.
%
% Inputs:
%   nt        : number of targets (integer)
%   zt(nt,1)  : target positions (complex)
%   p         : panel order (integer)
%   zsrc(p,1) : source nodes = s.x (complex)
%   nzsrc(p,1): source normals = s.nx (complex)
%   wzp(p,1)  : complex speed weights = s.wxp (complex)
%   a(1)      : panel start endpoint (complex scalar as 1-element array)
%   b(1)      : panel end   endpoint (complex scalar as 1-element array)
%   iside     : 1=exterior, 0=interior
%
% Outputs (p×nt, Fortran column-major):
%   As(p,nt)  : SLP weights (real double)
%   A(p,nt)   : DLP weights (complex double; take real for Laplace)
%   A1(p,nt)  : Re(Az) x-deriv (real)
%   A2(p,nt)  : -Im(Az) y-deriv (real)
%   A3(p,nt)  : Re(Azz) (real)
%   A4(p,nt)  : -Im(Azz) (real)
%
% Note: Fortran outputs are p×nt. MATLAB's SDspecialquad returns nt×p.
%   For comparison: As_matlab == As_fortran' (transpose).
%   For nt=1 (single target), both are length-p column vectors.
% --------------------------------------------------------------------------
