function [loop_tag] = gmsh_add_polygon_mex(n, pts_x, pts_y, mesh_size, loop_tag)
n         = double(n);
mesh_size = double(mesh_size);
loop_tag  = double(0);
mex_id_ = 'gmsh_add_polygon_(c i int64_t[x], c i double[x], c i double[x], c i double[x], c io int64_t[x])';
[loop_tag] = poisson2d_mex(mex_id_, n, pts_x, pts_y, mesh_size, loop_tag, 1, n, n, 1, 1);
end

% --------------------------------------------------------------------------
