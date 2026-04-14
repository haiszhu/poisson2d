classdef GmshGeo < handle
    % GmshGeo  MATLAB wrapper for Gmsh mesh generation.
    % Uses the mwrap-generated poisson2d_mex interface.
    %
    % Usage:
    %   geom = GmshGeo();
    %   geom.add_polygon(pts_Nx2, mesh_size);   % outer boundary
    %   geom.add_polygon(hole_Nx2, mesh_size);  % hole (reversed = CW)
    %   [nodes, tris] = geom.generate_mesh(2, 8);
    %   clear geom;  % calls delete -> gmsh_finalize_mex

    properties (Access = private)
        initialized = false;
    end

    methods
        function obj = GmshGeo()
            gmsh_init_mex();
            obj.initialized = true;
        end

        function loop_tag = add_polygon(obj, pts, mesh_size)
            % pts: Nx2 matrix of [x, y] coordinates
            % mesh_size: characteristic mesh length (scalar)
            n     = double(size(pts, 1));
            pts_x = double(pts(:, 1));
            pts_y = double(pts(:, 2));
            loop_tag = gmsh_add_polygon_mex(n, pts_x, pts_y, double(mesh_size));
        end

        function [nodes, tris] = generate_mesh(obj, dim, algorithm)
            % dim: must be 2 (only 2D supported)
            % algorithm: Gmsh meshing algorithm ID (default 8 = Frontal-Delaunay)
            if nargin < 2 || dim ~= 2
                error('GmshGeo:generate_mesh', 'Only 2D mesh supported.');
            end
            if nargin < 3, algorithm = 8; end
            [n_nodes, n_tris] = gmsh_generate_mex(double(algorithm));
            [nodes, tris]     = gmsh_fetch_mex(n_nodes, n_tris);
        end

        function delete(obj)
            if obj.initialized
                gmsh_finalize_mex();
                obj.initialized = false;
            end
        end
    end
end
