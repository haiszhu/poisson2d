function [z, zp, zpp, nz, w, wzp, awzp, zinter] = zinit_mex(npan, p, sinter, sinterdiff, T, W)
npan = double(npan);
p    = double(p);
np   = npan * p;
npi  = npan + 1;
z      = complex(zeros(np,  1));
zp     = complex(zeros(np,  1));
zpp    = complex(zeros(np,  1));
nz     = complex(zeros(np,  1));
wzp    = complex(zeros(np,  1));
zinter = complex(zeros(npi, 1));
w      = zeros(np,  1);
awzp   = zeros(np,  1);
mex_id_ = 'zinit_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c io dcomplex[x], c io dcomplex[x], c io dcomplex[x], c io dcomplex[x], c io double[x], c io dcomplex[x], c io double[x], c io dcomplex[x])';
[z, zp, zpp, nz, w, wzp, awzp, zinter] = poisson2d_mex(mex_id_, npan, p, sinter, sinterdiff, T, W, z, zp, zpp, nz, w, wzp, awzp, zinter, 1, 1, npi, npan, p, p, np, np, np, np, np, np, np, npi);
end

% --------------------------------------------------------------------------
