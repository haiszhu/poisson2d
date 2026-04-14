function [xq, yq, wq] = boxvolquad(BMmask, bounds, gridsize, p)
% BOXVOLQUAD  Tensor-product Chebyshev quadrature over active background boxes.
%
% [xq, yq, wq] = boxvolquad(BMmask, bounds, gridsize, p)
%
% Inputs:
%   BMmask   : N x N logical mask -- true for active boxes (from compositemesh)
%   bounds   : scalar, domain is [-bounds, bounds]^2
%   gridsize : scalar, box side length = 2*bounds/N
%   p        : Chebyshev order; p*p nodes per box
%
% Outputs:
%   xq  : (nboxes*p^2) x 1  x-coordinates
%   yq  : (nboxes*p^2) x 1  y-coordinates
%   wq  : (nboxes*p^2) x 1  weights  (sum over one box = gridsize^2)

[xc, wc, ~] = cheby(p);          % p nodes/weights on [-1,1]
[X, Y]   = meshgrid(xc, xc);     % p x p reference grids
[Wx, Wy] = meshgrid(wc, wc);
x_ref = X(:);                     % p^2 x 1
y_ref = Y(:);
w_ref = (Wx .* Wy)';
w_ref = w_ref(:);                 % p^2 x 1

h  = gridsize / 2;
[I, J] = find(BMmask);
nboxes = numel(I);
npp    = p * p;

xq = zeros(nboxes * npp, 1);
yq = zeros(nboxes * npp, 1);
wq = zeros(nboxes * npp, 1);

for k = 1:nboxes
    cx = -bounds + gridsize * (J(k) - 0.5);
    cy =  bounds - gridsize * (I(k) - 0.5);
    idx = (k-1)*npp + (1:npp);
    xq(idx) = cx + h * x_ref;
    yq(idx) = cy + h * y_ref;
    wq(idx) = h^2 * w_ref;
end
end
