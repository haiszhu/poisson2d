function [n_nodes, n_tris] = gmsh_generate_mex(algorithm, n_nodes, n_tris)
algorithm = double(algorithm);
n_nodes   = double(0);
n_tris    = double(0);
mex_id_ = 'gmsh_generate_(c i int64_t[x], c io int64_t[x], c io int64_t[x])';
[n_nodes, n_tris] = poisson2d_mex(mex_id_, algorithm, n_nodes, n_tris, 1, 1, 1);
end

% --------------------------------------------------------------------------
