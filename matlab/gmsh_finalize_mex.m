function [] = gmsh_finalize_mex()
mex_id_ = 'gmsh_finalize_()';
poisson2d_mex(mex_id_);
end

% --------------------------------------------------------------------------
% koorn_uvs_wts_mex: Vioreanu-Rokhlin nodes and weights on reference triangle
%   (0,0)-(1,0)-(0,1) for polynomial order norder.
%
% [uvs, wts] = koorn_uvs_wts_mex(norder)
%   norder       : integer order (0..20)
%   uvs(2,npols) : (u,v) node coordinates, npols=(norder+1)*(norder+2)/2
%   wts(npols,1) : integration weights, sum(wts)=0.5
% --------------------------------------------------------------------------
