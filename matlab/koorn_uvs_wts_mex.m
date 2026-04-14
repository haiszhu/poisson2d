function [uvs, wts] = koorn_uvs_wts_mex(norder, uvs, wts)
norder = double(norder);
npols  = (norder+1)*(norder+2)/2;
uvs    = zeros(2, npols);
wts    = zeros(npols, 1);
mex_id_ = 'koorn_uvs_wts_r64(c i int64_t[x], c io double[xx], c io double[x])';
[uvs, wts] = poisson2d_mex(mex_id_, norder, uvs, wts, 1, 2, npols, npols);
end
