function [xq, yq, wq] = trivolquad(nodes, tris, p)
% TRIVOLQUAD  Affine-mapped Vioreanu-Rokhlin quadrature over a triangle mesh.
%
% [xq, yq, wq] = trivolquad(nodes, tris, p)
%
% Inputs:
%   nodes  : Nv x 3  vertex coordinates (x, y, z) -- GmshGeo output
%   tris   : Nt x 3  triangle vertex indices, 1-based -- GmshGeo output
%   p      : polynomial order (integer, 0..20); npols=(p+1)*(p+2)/2 per tri
%
% Outputs:
%   xq     : (Nt*npols) x 1  x-coordinates of quadrature points
%   yq     : (Nt*npols) x 1  y-coordinates of quadrature points
%   wq     : (Nt*npols) x 1  quadrature weights (physical area element)
%
% Reference triangle: (0,0)-(1,0)-(0,1).  sum(wq over one triangle) = area.

[uvs, wts] = koorn_uvs_wts_mex(double(p));
npols = size(uvs, 2);   % (p+1)*(p+2)/2
u = uvs(1,:)';          % npols x 1
v = uvs(2,:)';          % npols x 1

ntris = size(tris, 1);
xq = zeros(ntris * npols, 1);
yq = zeros(ntris * npols, 1);
wq = zeros(ntris * npols, 1);

for k = 1:ntris
    x1 = nodes(tris(k,1), 1);  y1 = nodes(tris(k,1), 2);
    x2 = nodes(tris(k,2), 1);  y2 = nodes(tris(k,2), 2);
    x3 = nodes(tris(k,3), 1);  y3 = nodes(tris(k,3), 2);

    dx1 = x2 - x1;  dx2 = x3 - x1;
    dy1 = y2 - y1;  dy2 = y3 - y1;

    jacdet = abs(dx1*dy2 - dx2*dy1);   % constant for affine map

    idx = (k-1)*npols + (1:npols);
    xq(idx) = x1 + dx1*u + dx2*v;
    yq(idx) = y1 + dy1*u + dy2*v;
    wq(idx) = jacdet * wts;
end
end
