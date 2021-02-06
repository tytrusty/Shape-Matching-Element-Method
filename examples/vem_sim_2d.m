function vem_sim_2d
    % Example 2D VEM simulation.

    % Simulation parameters
    dt = 0.01;      	% timestep
    gravity = -10;     % gravity strength
    k_stability = 1e5;   % Stiffness for VEM stability term
    YM = 1e5; %in Pascals
    pr =  0.45;
    order = 1;          % (1 or 2) linear or quadratic deformation
    rho = 1e2 ;            % density
	d = 2;              % dimension (2 or 3)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    vol = 0.04;
    
    % The number of elements in the monomial basis.
    % -- example: Draft equation 1 is second order with k == 5
    k = basis_size(d, order);
    
    % Load a basic square mesh
    [V,F] = readOBJ('models/obj/plane.obj');
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
       E{i} = F(i,:)';
    end

    % Simpler cube example.
    x0  = [0 0; 0 2;  2 2; 2 0; ]';
    
    num_verts = 10;
    [X,Y] = meshgrid(linspace(0,2,num_verts), linspace(0,2,num_verts));
    V = [X(:) Y(:)];
    
    shapes_per_line = 4;
     
    % Generates line segment samples along a specified line.
    function s=sample_func(a,b)
        s1 = linspace(a(1),b(1),shapes_per_line+1);
        s2 = linspace(a(2),b(2),shapes_per_line+1);
        s = [];
        for ii=1:shapes_per_line
            s = [s [s1(ii:ii+1); s2(ii:ii+1)]];
        end
    end

    % Generate lines on each of the 4 lines of a box.
    xa = sample_func(x0(:,1),x0(:,2));
    xb = sample_func(x0(:,2),x0(:,3));
    xc = sample_func(x0(:,3),x0(:,4));
    xd = sample_func(x0(:,4),x0(:,1));
    x0 = [xa xb xc xd];
    E = reshape(1:size(x0,2), 2, []);
    F = E';
    E = num2cell(E, 1)';
        
    % min_I = find(x0(2,:) == max(x0(2,:))); % pin top side
    % min_I = find(x0(1,:) == min(x0(1,:))); % pin left side
    min_I = [1]; % Pinning the top left corner.

    % Initial deformed positions and velocities
    x = x0;
    v = zeros(size(x));
        
    % Setup figure
    fig=figure(1);
    clf;
    
    % Plot quadrature points
    X_plot = plot(V(:,1),V(:,2),'.');
    hold on;
    
    % Create plot of boundary edges.
    E_lines = cell(numel(E),1);
    cm=lines(numel(E));
    for i=1:numel(E)
        E_lines{i} = plot(x0(1,E{i}),x0(2,E{i}),'-','LineWidth',3, ...
            'Color', cm(i,:));
        hold on;
    end
    
    % Plot pinned vertices
    plot(x0(1,min_I), x0(2,min_I),'.','MarkerSize', 30,'Color','r');
    axis equal
    ylim([-4.5 2.5])
    
    % The number of quadrature points where at each point we
    % add stiffness matrix contributions.
    m = size(V,1);
    
    % Constraint matrix for boundary conditions.
    P = fixed_point_constraint_matrix(x0',sort(min_I)');
    
    % Lame parameters concatenated.
    [lambda, mu] = emu_to_lame(YM, pr);
    params = [mu * 0.5, lambda * 0.5];
    params = repmat(params,size(V,1),1);
    
    % Compute Shape weights0.8
    [w,w_I] = blending_weights_2d(V, x0', F, 1);
    [w0,w0_I] = blending_weights_2d(x0', x0', F, 1);
    
    % Generate centers of mass.
    [x0_coms, com_cluster, com_map] = generate_com(x0, E, w, numel(E));
    com_plt = plot(x0_coms(1,:),x0_coms(2,:), '.','Color','g','MarkerSize',20);
    hold on;
    
    % Build shape matching matrices
    L = compute_shape_matrices(x0, x0_coms, com_map, E, ...
        com_cluster, order, 'hierarchical');
    
    % Build weighted monomial bases
    [Y,Y_S,C_I] = vem_dx_dc(V', x0_coms, w, w_I, com_map, order, k);
    [Y0,Y0_S,C0_I] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, order, k);

    % Computing each gradient of deformation gradient with respect to
    % projection operator (c are polynomial coefficients)
    [dF_dc, dF_dc_S] = vem_dF_dc(V', x0_coms, w, w_I, com_map, order, k);

    % Gravity force vector.
    dg_dc = vem_ext_force([0 gravity]', repmat(rho .* vol, size(V,1),1), Y, Y_S);
    f_gravity = dt*P*(L' * dg_dc);

    % Compute mass matrices
%     ME = vem_error_matrix(Y0, L, w0_I, C0_I, d, k, numel(E));
%     M = vem_mass_matrix(Y, L, rho .* vol, w_I, C_I, d, k, numel(E));
    ME = vem_error_matrix_matlab(Y0, Y0_S, L, d);
    M = vem_mass_matrix_matlab(Y, Y_S, L, repmat(rho .* vol, size(V,1),1));
    M = (M + k_stability*ME); % sparse?

    ii=1;
    for t=0:dt:300
        b = [];
        for i=1:numel(E)
            b = [b x(:,E{i})];
        end
        b = b(:);
        
        % Solve for polynomial coefficients (projection operators).
        c = L * b;

        % force vector
        dV_dq = zeros(d*(k*numel(E) + size(x0_coms,2)),1); 
        
        % stiffness matrix
        K = zeros(d*(k*numel(E) + size(x0_coms,2)), d*(k*numel(E) + size(x0_coms,2)));
        
        % Computing force and stiffness contributions
        for i = 1:m
            % Deformation Gradient
            F = dF_dc{i} * dF_dc_S{i} * c;
            F = reshape(F,d,d);
            
            V(i,:) = Y{i} * Y_S{i} * c;
            
            % Force vector
            dV_dF = neohookean_dF(F, params(i,1), params(i,2));
            dV_dq = dV_dq +  dF_dc_S{i}' * dF_dc{i}' * dV_dF .* vol;
%             
            % Stiffness matrix contribution
            d2V_dF2 = neohookean_dF2(F, params(i,1), params(i,2));
            K = K - dF_dc_S{i}' * dF_dc{i}' * d2V_dF2 * dF_dc{i} * dF_dc_S{i} .* vol;
        end
        K = L' * K * L;
        dV_dq = L' * dV_dq;
        
        % Error correction force
        f_error = - 2 * ME * x(:);
        f_error = k_stability*(dt * P * f_error(:));
        
        % Force from potential energy.
        f_internal = -dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        lhs = P*(M - dt*dt*K)*P';
        rhs = P*M*v(:) + f_internal + f_gravity + f_error;
        qdot = lhs \ rhs;   % perform solve
        qdot = P'*qdot;     % unproject velocities
        v = reshape(qdot,2,[]);

        % Update position
        x = x + dt*v;

        % Update plots
        for i = 1:numel(E_lines)
            E_lines{i}.XData = x(1,E{i});
            E_lines{i}.YData = x(2,E{i});
        end
     
        X_plot.XData = V(:,1);
        X_plot.YData = V(:,2);
        drawnow
        
        if save_output
            fn=sprintf('output_png\\100k_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end