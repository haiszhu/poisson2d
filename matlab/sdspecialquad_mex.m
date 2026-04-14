function [As, Ad, A1, A2, A3, A4] = sdspecialquad_mex(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4)
nt    = double(nt);
p     = double(p);
iside = double(iside);
pnt   = p * nt;
if nargin < 10 || isempty(As),   As   = zeros(p, nt);           end
if nargin < 11 || isempty(Ad),   Ad   = complex(zeros(p, nt));  end
if nargin < 12 || isempty(A1),   A1   = zeros(p, nt);           end
if nargin < 13 || isempty(A2),   A2   = zeros(p, nt);           end
if nargin < 14 || isempty(A3),   A3   = zeros(p, nt);           end
if nargin < 15 || isempty(A4),   A4   = zeros(p, nt);           end
mex_id_ = 'sdspecialquad_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c io double[xx], c io dcomplex[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[As, Ad, A1, A2, A3, A4] = poisson2d_mex(mex_id_, nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4, 1, nt, 1, p, p, p, 1, 1, 1, p, nt, p, nt, p, nt, p, nt, p, nt, p, nt);
end

% ==========================================================================
% Gmsh mesh generation interface
% Backed by src/gmsh_wrap.cpp (extern "C" wrappers over Gmsh C++ API).
%
% Two-step generate+fetch pattern: gmsh_generate_mex returns counts n_nodes
% and n_tris; gmsh_fetch_mex fills pre-allocated double[n_nodes,3] arrays.
% GmshGeo.m calls these in sequence and hides the split from the user.
% ==========================================================================

% --------------------------------------------------------------------------
