% test_gmsh6.m - reproduce VIEsolver gmsh_utils/test_gmsh6.m in poisson2d
% Uses GmshGeo (mwrap poisson2d_mex interface) + compositemesh + quadr (local).
testdir = fileparts(mfilename('fullpath'));
addpath(fullfile(testdir, '../matlab'))

% 1. Define three curves (Outer Circle + Star 1 + Star 2)
sc1_val =  2.1 - 0.8*1i; r1_val = 1.0; a1_val = 0.25; f1_val = 3; rot1_val = 7*pi/5.0;
sc2_val = -2.0 + 1.3*1i; r2_val = 1.0; a2_val = 0.25; f2_val = 3; rot2_val = 2*pi/5.0;

curves{1}.Z   = @(t) 7 * exp(1i*t);
curves{1}.Zp  = @(t) 7i * exp(1i*t);
curves{1}.Zpp = @(t) -7 * exp(1i*t);
curves{1}.side = 'i';

% Star 1
curves{2}.Z   = @(t) sc1_val + (r1_val + r1_val*a1_val*cos(3*(t-rot1_val))) .* exp(1i*t);
curves{2}.Zp  = @(t) ((-3*r1_val*a1_val*sin(3*(t-rot1_val))) .* exp(1i*t)) + ...
                     ((r1_val + r1_val*a1_val*cos(3*(t-rot1_val))) .* (1i*exp(1i*t)));
curves{2}.Zpp = @(t) ((-9*r1_val*a1_val*cos(3*(t-rot1_val))) .* exp(1i*t)) + ...
                     ((-3*r1_val*a1_val*sin(3*(t-rot1_val))) .* (1i*exp(1i*t))) + ...
                     ((-3*r1_val*a1_val*sin(3*(t-rot1_val))) .* (1i*exp(1i*t))) + ...
                     ((r1_val + r1_val*a1_val*cos(3*(t-rot1_val))) .* (-exp(1i*t)));
curves{2}.side = 'e';

% Star 2
curves{3}.Z   = @(t) sc2_val + (r2_val + r2_val*a2_val*cos(3*(t-rot2_val))) .* exp(1i*t);
curves{3}.Zp  = @(t) ((-3*r2_val*a2_val*sin(3*(t-rot2_val))) .* exp(1i*t)) + ...
                     ((r2_val + r2_val*a2_val*cos(3*(t-rot2_val))) .* (1i*exp(1i*t)));
curves{3}.Zpp = @(t) ((-9*r2_val*a2_val*cos(3*(t-rot2_val))) .* exp(1i*t)) + ...
                     ((-3*r2_val*a2_val*sin(3*(t-rot2_val))) .* (1i*exp(1i*t))) + ...
                     ((-3*r2_val*a2_val*sin(3*(t-rot2_val))) .* (1i*exp(1i*t))) + ...
                     ((r2_val + r2_val*a2_val*cos(3*(t-rot2_val))) .* (-exp(1i*t)));
curves{3}.side = 'e';

%
[s2 N] = quadr(curves{2}, 256, 'p', 'g');
[s3 N] = quadr(curves{3}, 256, 'p', 'g');
A2 = 0.5 * imag(sum(conj(s2.x - mean(s2.x)) .* s2.xp .* s2.w));
A3 = 0.5 * imag(sum(conj(s3.x - mean(s3.x)) .* s3.xp .* s3.w));

% 2. Kress nodes
NBS = 320;
tt = linspace(0, 2*pi, NBS+1); tt(end) = [];
for ell = 1:3
    kress_nodes{ell} = curves{ell}.Z(tt);
end

% 3. compositemesh
bounds = 13; Nmesh = 64;
fprintf('Generating all polygons*.mat files...\n');
[BMmask, polygons0, polygonsi] = compositemesh(curves, bounds, Nmesh, kress_nodes);
fprintf('Done.\n');

% 4. Triangulation, Plotting, and Quadrature Assembly
figure(100); clf; hold on;
colors = {'cyan', 'yellow', 'magenta'};

gridsize = 2 * bounds / Nmesh;
p = 10;   % quadrature order for all three region types

xq_tris = []; yq_tris = []; wq_tris = [];   % straight triangles (Gmsh)
xq_ctri = []; yq_ctri = []; wq_ctri = [];   % curved triangles (Gordon-Hall VR)
xq_box  = []; yq_box  = []; wq_box  = [];   % background boxes

% Outer boundary region
j = 0;
fprintf('--- Generating Triangle Mesh for Region %d ---\n', j);
geom = GmshGeo();
geom.add_polygon(polygons0.pOarr, polygons0.gridsize*2);
geom.add_polygon(polygons0.pIarr, polygons0.gridsize*2);
[nodes, tris] = geom.generate_mesh(2, 8);
clear geom;
patch('Faces', tris, 'Vertices', nodes(:,1:2), ...
      'FaceColor', colors{j+1}, 'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 0.6);
plot(polygons0.pOarr(:,1), polygons0.pOarr(:,2), 'k--', 'LineWidth', 1);
plot(polygons0.pIarr(:,1), polygons0.pIarr(:,2), 'k--', 'LineWidth', 1);
[xq, yq, wq] = trivolquad(nodes, tris, p);
xq_tris = [xq_tris; xq]; yq_tris = [yq_tris; yq]; wq_tris = [wq_tris; wq];
% curved triangles for outer boundary
xq_ctri = [xq_ctri; polygons0.ctrisXq(:)];
yq_ctri = [yq_ctri; polygons0.ctrisYq(:)];
wq_ctri = [wq_ctri; polygons0.ctrisWq(:)];

% Interior islands
for j = 1:2
  fprintf('--- Generating Triangle Mesh for Region %d ---\n', j);
  geom = GmshGeo();
  geom.add_polygon(polygonsi{j}.pOarr, polygonsi{j}.gridsize*2);
  geom.add_polygon(polygonsi{j}.pIarr, polygonsi{j}.gridsize*2);
  [nodes, tris] = geom.generate_mesh(2, 8);
  clear geom;
  patch('Faces', tris, 'Vertices', nodes(:,1:2), ...
        'FaceColor', colors{j+1}, 'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 0.6);
  plot(polygonsi{j}.pOarr(:,1), polygonsi{j}.pOarr(:,2), 'k--', 'LineWidth', 1);
  plot(polygonsi{j}.pIarr(:,1), polygonsi{j}.pIarr(:,2), 'k--', 'LineWidth', 1);
  [xq, yq, wq] = trivolquad(nodes, tris, p);
  xq_tris = [xq_tris; xq]; yq_tris = [yq_tris; yq]; wq_tris = [wq_tris; wq];
  % curved triangles for island j
  xq_ctri = [xq_ctri; polygonsi{j}.ctrisXq(:)];
  yq_ctri = [yq_ctri; polygonsi{j}.ctrisYq(:)];
  wq_ctri = [wq_ctri; polygonsi{j}.ctrisWq(:)];
end

% Background boxes
[xq_box, yq_box, wq_box] = boxvolquad(BMmask, bounds, gridsize, p);

% Draw regular grid boxes
[I, J] = find(BMmask);
for k = 1:length(I)
  cx = -bounds + gridsize*(J(k)-0.5);
  cy =  bounds - gridsize*(I(k)-0.5);
  bx = cx + [-gridsize/2, gridsize/2, gridsize/2, -gridsize/2, -gridsize/2];
  by = cy + [ gridsize/2,  gridsize/2, -gridsize/2, -gridsize/2,  gridsize/2];
  plot(bx, by, 'k--', 'LineWidth', 0.1);
end

fprintf('Quadrature point counts:\n');
fprintf('  curved triangles : %d pts\n', numel(wq_ctri));
fprintf('  straight triangles: %d pts\n', numel(wq_tris));
fprintf('  background boxes : %d pts\n', numel(wq_box));
fprintf('Area check (integral of 1):\n');
fprintf('  curved tri  sum(wq) = %.6f\n', sum(wq_ctri));
fprintf('  straight tri sum(wq) = %.6f\n', sum(wq_tris));
fprintf('  box         sum(wq) = %.6f\n', sum(wq_box));
fprintf('  total               = %.6f\n', sum(wq_ctri)+sum(wq_tris)+sum(wq_box));

% Original smooth curves
t_plot = linspace(0, 2*pi, 1000);
for ell = 1:3
  Z_val = curves{ell}.Z(t_plot);
  plot(real(Z_val), imag(Z_val), 'r-', 'LineWidth', 2);
end

axis equal; grid on;
title('poisson2d GmshGeo: CompositeMesh Triangulation');
xlabel('x'); ylabel('y');
saveas(gcf, 'test_compositemesh.png');
fprintf('Saved test_gmsh6.png\n');

figure(101),clf,
plot(xq_ctri,yq_ctri,'.'); axis equal, hold on
plot(xq_tris,yq_tris,'.');
plot(xq_box,yq_box,'.');

area_exact = 49*pi - A2 - A3;
area_quad  = sum(wq_ctri) + sum(wq_tris) + sum(wq_box);
fprintf('Area verification:\n');
fprintf('  exact  (49*pi - A2 - A3) = %.15f\n', area_exact);
fprintf('  quad   (ctri+tri+box)    = %.15f\n', area_quad);
fprintf('  error                    = %.6e\n',  area_exact - area_quad);