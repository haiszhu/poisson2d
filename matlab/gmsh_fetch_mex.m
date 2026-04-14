function [nodes, tris] = gmsh_fetch_mex(n_nodes, n_tris, nodes, tris)
n_nodes = double(n_nodes);
n_tris  = double(n_tris);
nodes   = zeros(n_nodes, 3);
tris    = zeros(n_tris,  3);
mex_id_ = 'gmsh_fetch_(c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx])';
[nodes, tris] = poisson2d_mex(mex_id_, n_nodes, n_tris, nodes, tris, 1, 1, n_nodes, 3, n_tris, 3);
end

% --------------------------------------------------------------------------
