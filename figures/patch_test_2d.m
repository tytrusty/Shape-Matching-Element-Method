function patch_test_2d
% Simulation parameters
    order = 1;          % (1 or 2) linear or quadratic deformation
	d = 2;              % dimension (2 or 3)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    
    figure(1);
    clf; axis equal;
    % The number of elements in the monomial basis.
    % -- example: Draft equation 1 is second order with k == 5
    k = basis_size(d, order);
    
    % Load a basic square mesh
    [V,F,UV,TF] = readOBJ('models/obj/plane.obj');
    F0=F;
    V0=V;
    z = V(:,3);
    V = V(:,1:2);
    
    % Get boundary edges. The points on each edge will be our nodal
    % values used in computing the "projection" operator in VEM.
    % Each edge here represents a single "shape" in the shape matching
    % context.
    F = boundary_faces(F);
        
    % Get undeformed boundary vertices. 
    V_bnd = unique(F(:));
    x0 = V(V_bnd,:)';

    % Remap edge vertex IDs to IDs into 'x0'
    new_Idx = find(V_bnd);
    [toreplace, bywhat] = ismember(F, V_bnd);
    F(toreplace) = new_Idx(bywhat(toreplace));
    
    % Todo -- only use cells for all edge stuff
    E = cell(size(F,1),1);
    for i =1:size(F,1)
       E{i} = F(i,:);
    end
    
    % Plot quadrature points
    X_plot = plot(V(:,1),V(:,2),'.','Color','g','MarkerSize',10);
    hold on;
    
    % Create plot of boundary edges.
    E_lines = cell(numel(E),1);
    cm=lines(numel(E));
    cm=copper(numel(E));
    for i=1:numel(E)
        E_lines{i} = plot(x0(1,E{i}),x0(2,E{i}),'-','LineWidth',2, ...
            'Color', cm(i,:));
        hold on;
    end
    
    % Undeformed Center of mass
    % -- Draft equation (4)
    x0_com = mean(x0,2);
    com_map = ones(numel(E),1);
    com_cluster = {}; % TODO: verify this can be empty for global
    
    % Build shape matching matrices
    L = compute_shape_matrices(x0, x0_com, com_map, E, ...
        com_cluster, order, 'global');
    
    x = x0;
    % Linear transform
%     x(1,:) = 0.8 * x(1,:);
%     x(2,:) = (3/4) * x(2,:);
%         x(2,:) = 1.5 * x(2,:);

    % Rigid transform
    theta = 0;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    R(1,2) = 0.5;
%     R(2,1) = 2;

    x = R * x;
    
    E_lines2 = cell(numel(E),1);
    cm=copper(numel(E));
    for i=1:numel(E)
        E_lines2{i} = plot(x(1,E{i}),x(2,E{i}),'-','LineWidth',2, ...
            'Color', cm(i,:));
        hold on;
    end
    
    % Shape matching
    b = [];
    for i=1:numel(E)
        b = [b (x(:,E{i}))];
    end
    b = b(:);
    c = L * b;
    x0=x0'
    function [a,w] = weights(X)
        dE = (V0(F(:,2),:) - V0(F(:,1),:));
        norm2_dE = sum(dE.*dE, 2);

        dVx = X(:,1)' - V0((F(:,1)),1);
        dVy = X(:,2)' - V0((F(:,1)),2);

        dEdotdV = (dE(:,1).*dVx +  dE(:,2).*dVy)./norm2_dE;

        %bound constraints
        dEdotdV = min(max(dEdotdV, 0),1);

        a = zeros(size(X,1), size(F,1));
        w = zeros(size(X,1), size(F,1));

        %%% THIS BIT DOES THE RAYTRACING WEIGHTS %%%%%
        for jj=1:size(F,1)

            %nearest point on i^th edge of j^th point (column)
            VE_X = (dEdotdV(jj,:).*dE(jj,1) + x0(F(jj,1),1))';
            VE_Y = (dEdotdV(jj,:).*dE(jj,2) + x0(F(jj,1),2))';

            %build ray from nearest point to i^th point 
            R_X = X(:,1) - VE_X;
            R_Y = X(:,2) - VE_Y;

            dE_X = dE(:,1)';
            dE_Y = dE(:,2)';

            %intersect every edge with very edge 
            det_A = 1./(-R_Y.*dE_X + R_X.*dE_Y);

            %coefficient along the plane
            B_X = X(:,1) - x0(F(:,1),1)';
            B_Y = X(:,2) - x0(F(:,1),2)';

            first_row = -R_Y.*B_X + R_X.*B_Y;
            second_row = dE_X.*B_Y - dE_Y.*B_X;

            %intersection parameter
            beta = first_row.*det_A;
            gamma = second_row.*det_A;

            beta(gamma < 0) = inf;
            gamma(gamma < 0) = inf;

            %beta(:,jj) = inf(size(V,1),1);
            %gamma(:,jj) = inf(size(V,1),1);

            gamma(beta > 1) = inf;
            gamma(beta < 0) = inf;

            beta(beta > 1) = inf;
            beta(beta < 0) = inf;

            %deal with special cases
            beta(isinf(det_A)) = inf;
            beta(and(abs(R_X) ==0, abs(R_Y) == 0)) = 0;

            gamma(isinf(det_A)) = inf;
            gamma(and(abs(R_X) == 0, abs(R_Y) == 0)) = 0;

            %reconstruct the points I hit
            [min_gamma, edge_id] = min(gamma, [], 2);

            hit_points = diag(beta(sub2ind(size(beta), (1:size(X,1))', edge_id)))*dE(edge_id,1:2)  + V0(F(edge_id,1),1:2);

            dt = sqrt((hit_points(:,1) - VE_X).^2 + (hit_points(:,2) - VE_Y).^2);
            dx = sqrt((X(:,1) - VE_X).^2 + (X(:,2) - VE_Y).^2);
    %         w(:,jj) = max(1.0 - dx./0.8,0);
            w(:,jj) = max(1.0 - dx./min(0.5,dt),0);
        end

        for ii = 1:size(a,1)

            W = diag(1./w(ii,:));
            W(isinf(W)) = 1e8;

            a(ii,:) = quadprog(W, zeros(size(F,1),1), [], [],...
                                ones(1,size(F,1)), 1,...
                               zeros(size(F,1),1), ones(size(F,1),1));

        end
    %     w = a;
    end
%     w = compute_projected_weights(x0, E, V');
    w = weights(V);
    w_I = cell(size(V,1),1);
    for i=1:size(V,1)
        w_I{i} = 1:numel(E);
    end
    
    % Build Monomial bases for all quadrature points
    [Y,~] = vem_dx_dc(V', x0_com, w, w_I, com_map, order, k);
    
    for i=1:size(V,1)
        V(i,:) = (Y{i} * c)';
    end
    
    V_plot2 = plot(V(:,1),V(:,2),'.','Color','r','MarkerSize',10);
    hold on;
    V = [V z];
    writeOBJ('plane_output.obj', V,F0,UV,TF);
    axis equal;
end