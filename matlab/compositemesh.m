function [BMmask,polygons0,polygonsi] = compositemesh(curves, bounds, N, kress_nodes)
    % MATLAB implementation of CompositeMesh.__init__ logic
    % Replicates the generation of polygons*.mat files
    % 
    % curves: cell array {Z, Zp, Zpp, 'i'/'e'}
    % bounds: scalar (defines [-bounds, bounds]^2)
    % N: grid resolution
    % kress_nodes: cell array of collocation points for each curve
    
    nbdry = length(curves);
    gridsize = 2 * bounds / N;
    order = 16;
    [xg,wg] = gauss(order);
    
    % Background grid points (centroids)
    [XX, YY] = meshgrid(linspace(-bounds + gridsize/2, bounds - gridsize/2, N), ...
                        linspace(bounds - gridsize/2, -bounds + gridsize/2, N));
    sclist = XX + 1j*YY;
    
    % Background grid boxes for polyshape
    % We'll use these for the unary_union equivalent
    
    BMmask = true(N, N);
    meshagg = cell(nbdry, 1);
    
    for ell = 1:nbdry
        curve = curves{ell};
        side = curve.side;
        
        nt = 10*numel(kress_nodes{ell}); % very rough
        % tt = linspace(0.05, 2*pi + 0.05, nt+1); % what is this... in python
        tt = linspace(0, 2*pi, nt+1);
        s = quadr(struct('x',curve.Z(tt(1:end-1))));
        ws = s.ws;
        dl = sum(s.ws) / ceil(sum(s.ws) / gridsize);
        Z = s.Z;
        Zp = s.Zp;
        
        % normal shooting points (norms_t)
        norms_t_ind = zeros(1, nt + 1); 
        norms_t_ind(1) = 1; 
        nElem = 1;
        acc = 0;
        total_len = sum(ws);
        curr_len = 0;
        for j = 1:nt-1
          curr_len = curr_len + ws(j);
          acc = acc + ws(j);
          if acc >= dl
            if (total_len - curr_len) >= 0.8 * dl % add if remaining space to the end is at least 80% of dl
              nElem = nElem + 1;
              norms_t_ind(nElem) = j + 1; 
              acc = 0;
            else % < 0.8 dl, go back and take average, worse case would be ~50% 
              nElem = nElem + 1;
              norms_t_ind(nElem) = round((nt+1 + norms_t_ind(nElem-1))/2);
              break
            end
          end
        end
        %
        nElem = nElem + 1;
        norms_t_ind(nElem) = nt + 1;
        norms_t_ind = norms_t_ind(1:nElem); % truncate
        norms_t = tt(norms_t_ind);
        nElem = nElem-1;
        
        % For calculating buffer zones
        norms_x = Z(norms_t);
        midpts_t = (norms_t(1:end-1) + norms_t(2:end)) / 2;
        nodePts = zeros(2*nElem,1);
        for j = 1:nElem
          % 
          dist = abs(norms_x(j) - norms_x(j+1));
          xp_mid = Zp(midpts_t(j));
          nxloc = -1j * xp_mid / abs(xp_mid);
          %
          if ell == 1 % Outer boundary (Python ell=0)
            nodePts(2*j-1) = norms_x(j);
            nodePts(2*j)   = Z(midpts_t(j)) - min([dist, gridsize, dl]) * nxloc;
          else
            nodePts(2*j-1) = norms_x(j);
            nodePts(2*j) = Z(midpts_t(j)) + min([dist, gridsize, dl]) * nxloc;
          end
        end
        
        ctrisNodePts = zeros(3,nElem);
        ctrisSidePts = zeros(order,3,nElem);
        % For curved triangle: nodes and GL-sampled side points
        %   side 1: straight v1 -> v2  (curve pt -> shot pt)
        %   side 2: straight v2 -> v3  (shot pt -> next curve pt)
        %   side 3: curved  v3 -> v1   (arc of Z, from norms_t(j+1) to norms_t(j))
        if ell == 1
          for j = 1:nElem-1
            ctrisNodePts(1,j) = nodePts(2*j-1);     % curve pt at norms_t(j)
            ctrisNodePts(2,j) = nodePts(2*j);       % shot pt
            ctrisNodePts(3,j) = nodePts(2*j+1);     % curve pt at norms_t(j+1)
            ctrisSidePts(:,1,j) = nodePts(2*j-1) + (xg+1)/2 * (nodePts(2*j)   - nodePts(2*j-1));
            ctrisSidePts(:,2,j) = nodePts(2*j)   + (xg+1)/2 * (nodePts(2*j+1) - nodePts(2*j));
            t_GL = norms_t(j+1) + (xg+1)/2 * (norms_t(j) - norms_t(j+1));
            ctrisSidePts(:,3,j) = Z(t_GL);
          end
          % Last triangle: v3 wraps to first curve point
          j = nElem;
          ctrisNodePts(1,j) = nodePts(2*j-1);
          ctrisNodePts(2,j) = nodePts(2*j);
          ctrisNodePts(3,j) = nodePts(1);
          ctrisSidePts(:,1,j) = nodePts(2*j-1) + (xg+1)/2 * (nodePts(2*j) - nodePts(2*j-1));
          ctrisSidePts(:,2,j) = nodePts(2*j)   + (xg+1)/2 * (nodePts(1)   - nodePts(2*j));
          t_GL = norms_t(j+1) + (xg+1)/2 * (norms_t(j) - norms_t(j+1));
          ctrisSidePts(:,3,j) = Z(t_GL);
        else
          for j = 1:nElem-1
            ctrisNodePts(1,j) = nodePts(2*j-1);
            ctrisNodePts(2,j) = nodePts(2*j);
            ctrisNodePts(3,j) = nodePts(2*j+1);
            ctrisSidePts(:,1,j) = nodePts(2*j-1) + (xg+1)/2 * (nodePts(2*j)   - nodePts(2*j-1));
            ctrisSidePts(:,2,j) = nodePts(2*j)   + (xg+1)/2 * (nodePts(2*j+1) - nodePts(2*j));
            t_GL = norms_t(j) + (xg+1)/2 * (norms_t(j+1) - norms_t(j));
            ctrisSidePts(:,3,j) = Z(t_GL(end:-1:1));
          end
          j = nElem;
          ctrisNodePts(1,j) = nodePts(2*j-1);
          ctrisNodePts(2,j) = nodePts(2*j);
          ctrisNodePts(3,j) = nodePts(1);
          ctrisSidePts(:,1,j) = nodePts(2*j-1) + (xg+1)/2 * (nodePts(2*j) - nodePts(2*j-1));
          ctrisSidePts(:,2,j) = nodePts(2*j)   + (xg+1)/2 * (nodePts(1)   - nodePts(2*j));
          t_GL = norms_t(j) + (xg+1)/2 * (norms_t(j+1) - norms_t(j));
          ctrisSidePts(:,3,j) = Z(t_GL(end:-1:1));
        end
        
        % --- Gordon-Hall mapped Vioreanu-Rokhlin volume quadrature (all ell) ---
        % ctrisNodePts convention: row1=curve pt at norms_t(j) [xi=0],
        %                          row2=shot pt, row3=curve pt at norms_t(j+1) [xi=1]
        % Gordon-Hall chi(xi)=Z(norms_t(j)+xi*dt), same formula for all ell.
        % ell=1: outer boundary, triangles CCW -> jac>0
        % ell>1: inner holes, orientation may be CW -> use abs(jac)
        p_vr = 10;
        [uvs_vr, wts_vr] = koorn_uvs_wts_mex(p_vr);
        npols = size(uvs_vr, 2);            % (p_vr+1)*(p_vr+2)/2 = 66
        xi_k  = uvs_vr(1,:);               % 1 x npols
        eta_k = uvs_vr(2,:);               % 1 x npols
        ctrisXq = zeros(npols, nElem);
        ctrisYq = zeros(npols, nElem);
        ctrisWq = zeros(npols, nElem);
        for j = 1:nElem
            % Vertices: x1,y1 = curve pt at xi=0; x2,y2 = curve pt at xi=1; x3,y3 = shot pt
            x1 = real(ctrisNodePts(1,j));  y1 = imag(ctrisNodePts(1,j));
            x2 = real(ctrisNodePts(3,j));  y2 = imag(ctrisNodePts(3,j));
            x3 = real(ctrisNodePts(2,j));  y3 = imag(ctrisNodePts(2,j));
            % Curved edge: Z parametrized forward xi=0->1 = t_j->t_{j+1}
            dt    = norms_t(j+1) - norms_t(j);
            t_k   = norms_t(j) + xi_k * dt;            % 1 x npols
            chi   = real(Z(t_k));                        % 1 x npols
            pi_   = imag(Z(t_k));                        % 1 x npols
            dZdt  = Zp(t_k) * dt;                       % 1 x npols, d(Z)/d(xi)
            chi_d = real(dZdt);  pi_d = imag(dZdt);
            % Deviation of curved edge from affine edge x1->x2
            f_c = chi - ((1-xi_k)*x1 + xi_k*x2);       % 1 x npols
            f_p = pi_ - ((1-xi_k)*y1 + xi_k*y2);       % 1 x npols
            s   = (1-xi_k-eta_k) ./ (1-xi_k);          % blend factor, 1 x npols
            % Gordon-Hall map to physical coordinates
            xq = (1-xi_k-eta_k)*x1 + xi_k*x2 + eta_k*x3 + f_c.*s;
            yq = (1-xi_k-eta_k)*y1 + xi_k*y2 + eta_k*y3 + f_p.*s;
            % Jacobian: ds/dxi = -eta/(1-xi)^2, ds/deta = -1/(1-xi)
            dxdxi  = (x2-x1) + (chi_d-(x2-x1)).*s - f_c.*eta_k./(1-xi_k).^2;
            dxdeta = (x3-x1) - f_c./(1-xi_k);
            dydxi  = (y2-y1) + (pi_d-(y2-y1)).*s - f_p.*eta_k./(1-xi_k).^2;
            dydeta = (y3-y1) - f_p./(1-xi_k);
            jac    = dxdxi.*dydeta - dxdeta.*dydxi;
            ctrisXq(:,j) = xq(:);
            ctrisYq(:,j) = yq(:);
            if ell == 1
                ctrisWq(:,j) = wts_vr .* jac(:);        % CCW -> jac > 0
            else
                ctrisWq(:,j) = wts_vr .* abs(jac(:));   % holes may be CW
            end
        end

        % Identify grid boxes involved in the buffer zone as false
        cmask = true(N, N);
        corners = [ -1 -1;...
                     1 -1;...
                     1  1;...
                    -1  1 ] * (gridsize/2);
        for i = 1:N
          for k = 1:N
            pt = sclist(i,k);
            corner_pts = pt + corners(:,1) + 1j*corners(:,2);
            min_dist = min(min(abs(corner_pts - nodePts(:).')));
            if (min_dist < 0.8 * gridsize) % if close to zig-zag nodePts
              cmask(i,k) = false;
            end
          end
        end
        
        % BMmask is what we want
        nodeX = real(nodePts);
        nodeY = imag(nodePts);
        if ell == 1 % interior, side
            % box inside zig-zag?
            [IN, ~] = inpolygon(real(sclist), imag(sclist), nodeX, nodeY);
            cmask = (~cmask) & IN;
            bmask = IN & (~cmask);
            [I, J] = find(bmask); % for k=1:numel(I), plot(XX(I(k),J(k)), YY(I(k),J(k)), '*'); end
            polygons0.bmask = bmask;

            % Python flips if NOT CCW for interior (wait, let's check)
            % s1 has sides='i', Python ell=0 logic
            nodeOx = nodeX; nodeOy = nodeY;
            if nodes_are_ccw(nodeOx, nodeOy)
              nodeOx = fliplr(nodeOx);
              nodeOy = fliplr(nodeOy);
            end
            
            merged = union_boxes(I, J, bounds, gridsize, N);
            [nodeIx, nodeIy] = boundary(merged);
            nodeIx = nodeIx(1:end-1).'; nodeIy = nodeIy(1:end-1).'; % Col vector to row
            if nodes_are_ccw(nodeIx, nodeIy)
                nodeIx = fliplr(nodeIx);
                nodeIy = fliplr(nodeIy);
            end
            
            % Update global grid mask (keep points INSIDE outer boundary)
            [IN_grid, ~] = inpolygon(real(sclist), imag(sclist), [nodeIx, nodeIx(1)], [nodeIy, nodeIy(1)]);
            BMmask = BMmask & IN_grid;
        else
            % box outside zig-zag?
            [IN, ~] = inpolygon(real(sclist), imag(sclist), nodeX, nodeY);
            cmask = (~cmask) & (~IN);
            bmask = (~IN) & (~cmask);
            [I, J] = find(bmask); % for k=1:numel(I), plot(XX(I(k),J(k)), YY(I(k),J(k)), '*'); end
            polygonsi{ell-1}.bmask = bmask;

            nodeIx = nodeX; nodeIy = nodeY;
            if nodes_are_ccw(nodeIx, nodeIy)
                nodeIx = fliplr(nodeIx);
                nodeIy = fliplr(nodeIy);
            end
            
            merged = union_boxes(I, J, bounds, gridsize, N);
            [nodeOx, nodeOy] = boundary(merged, 2); 
            nodeOx = nodeOx(1:end-1).'; nodeOy = nodeOy(1:end-1).';
            if nodes_are_ccw(nodeOx, nodeOy)
                nodeOx = fliplr(nodeOx);
                nodeOy = fliplr(nodeOy);
            end
            
            % Update global grid mask (remove points INSIDE holes)
            [IN_grid, ~] = inpolygon(real(sclist), imag(sclist), [nodeOx, nodeOx(1)], [nodeOy, nodeOy(1)]);
            BMmask = BMmask & (~IN_grid);
        end
        
        meshagg{ell}.nodeIx       = nodeIx;
        meshagg{ell}.nodeIy       = nodeIy;
        meshagg{ell}.nodeOx       = nodeOx;
        meshagg{ell}.nodeOy       = nodeOy;
        meshagg{ell}.bmask        = bmask;
        meshagg{ell}.norms_t      = norms_t;
        meshagg{ell}.norms_t_ind  = norms_t_ind;
        meshagg{ell}.ctrisNodePts = ctrisNodePts;
        meshagg{ell}.ctrisSidePts = ctrisSidePts;
        meshagg{ell}.ctrisXq      = ctrisXq;
        meshagg{ell}.ctrisYq      = ctrisYq;
        meshagg{ell}.ctrisWq      = ctrisWq;
    end
    
    % Now perform the actual meshing and save mat files
    % polygons0.mat: outer circle gap
    mesh0 = meshagg{1};
    pOarr = [mesh0.nodeOx(:), mesh0.nodeOy(:)];
    pIarr = [mesh0.nodeIx(:), mesh0.nodeIy(:)];
    norms_t = mesh0.norms_t; 
    norms_t_ind = mesh0.norms_t_ind;
    save('polygons0.mat', 'pOarr', 'pIarr', 'gridsize', 'norms_t', 'norms_t_ind', 'BMmask');
    meshagg(1) = []; % Remove processed
    polygons0.pOarr        = pOarr;
    polygons0.pIarr        = pIarr;
    polygons0.gridsize     = gridsize;
    polygons0.norms_t      = norms_t;
    polygons0.norms_t_ind  = norms_t_ind;
    polygons0.ctrisNodePts = mesh0.ctrisNodePts;
    polygons0.ctrisSidePts = mesh0.ctrisSidePts;
    polygons0.ctrisXq      = mesh0.ctrisXq;
    polygons0.ctrisYq      = mesh0.ctrisYq;
    polygons0.ctrisWq      = mesh0.ctrisWq;
    
    island_count = 1;
    while ~isempty(meshagg)
        mesh_master = meshagg{1};
        meshagg(1) = [];
        
        % For islands
        pOarr = [mesh_master.nodeOx(:), mesh_master.nodeOy(:)];
        pIarr = [mesh_master.nodeIx(:), mesh_master.nodeIy(:)];
        norms_t = mesh_master.norms_t;
        norms_t_ind = mesh_master.norms_t_ind;
        save(sprintf('polygons%d.mat', island_count), 'pOarr', 'pIarr', 'gridsize', 'norms_t', 'norms_t_ind');
        
        polygonsi{island_count}.pOarr        = pOarr;
        polygonsi{island_count}.pIarr        = pIarr;
        polygonsi{island_count}.gridsize     = gridsize;
        polygonsi{island_count}.norms_t      = norms_t;
        polygonsi{island_count}.norms_t_ind  = norms_t_ind;
        polygonsi{island_count}.ctrisNodePts = mesh_master.ctrisNodePts;
        polygonsi{island_count}.ctrisSidePts = mesh_master.ctrisSidePts;
        polygonsi{island_count}.ctrisXq      = mesh_master.ctrisXq;
        polygonsi{island_count}.ctrisYq      = mesh_master.ctrisYq;
        polygonsi{island_count}.ctrisWq      = mesh_master.ctrisWq;

        island_count = island_count + 1;
    end
end

function tf = nodes_are_ccw(x, y)
    % Signed area using shoelace formula
    % Area > 0 for CCW
    if isempty(x), tf = false; return; end
    nx = x(:); ny = y(:);
    area = sum(nx .* [ny(2:end); ny(1)] - [nx(2:end); nx(1)] .* ny);
    tf = area > 0;
end

function merged = union_boxes(I, J, bounds, gridsize, N)
    % Helper to union grid boxes using polyshape
    if isempty(I)
        merged = polyshape();
        return;
    end
    
    ps = repmat(polyshape(), 1, length(I));
    for k = 1:length(I)
        cx = -bounds + gridsize*(J(k)-0.5);
        cy =  bounds - gridsize*(I(k)-0.5);
        px = cx + [-gridsize/2, gridsize/2, gridsize/2, -gridsize/2];
        py = cy + [ gridsize/2,  gridsize/2, -gridsize/2, -gridsize/2];
        ps(k) = polyshape(px, py);
    end
    merged = union(ps);
end
