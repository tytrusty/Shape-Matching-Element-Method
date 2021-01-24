function patch_test_2d
% Simulation parameters
    order = 1;          % (1 or 2) linear or quadratic deformation
	d = 2;              % dimension (2 or 3)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    
    figure(1);
    clf;
    % The number of elements in the monomial basis.
    % -- example: Draft equation 1 is second order with k == 5
    k = basis_size(d, order);
    
    % Load a basic square mesh
    [V,F,UV,TF] = readOBJ('models/plane3.obj');
    F0=F;
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
    x(1,:) = 2 * x(1,:);
%     x(2,:) = (3/4) * x(2,:);
    
    % Rigid transform
%     theta = 45;
%     R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%     x = R * x;
    
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
    
    w = compute_projected_weights(x0, E, V');
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
    
end