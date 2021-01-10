function vem_sim_2d
    % Example 2D VEM simulation.

    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 170;   	% Lame parameter 1
    D = 0.5 * 1500;   	% Lame parameter 2
    gravity = -100;     % gravity strength
    k_error = 100000;   % Stiffness for VEM stability term
    order = 1;          % (1 or 2) linear or quadratic deformation
    rho = 1;            % density
	d = 2;              % dimension (2 or 3)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    
    % The number of elements in the monomial basis.
    % -- example: Draft equation 1 is second order with k == 5
    k = basis_size(d, order);
    
    % Load a basic square mesh
    [V,F] = readOBJ('models/plane.obj');
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

    % Simpler cube example.
    x0  = [0 0; 0 2;  2 2; 2 0; ]';
    
    num_verts = 10;
    [X,Y] = meshgrid(linspace(0,2,num_verts), linspace(0,2,num_verts));
    V = [X(:) Y(:)];
    
    shapes_per_line = 1;
     
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
    E = num2cell(E', 2);

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
    
    % Gravity force vector.
  	f_gravity = repmat([0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    % -- Draft equation (4)
    x0_com = mean(x0,2);
    
    % Build shape matching matrices
    %
    % In Shape matching paper, the projection operator is defined as
    % A = A_pq * A_qq (Shape matching equation 7)
    % The 'p' part is the only part that changes during simulation,
    % so we instead write A = PB where P = [p_1 p_2 ... p_n].
    %
    % Note: each shape (edge in this case) has its own 'B' matrix
    %       because each shape has its own projection operator.
    L = compute_shape_matrices(x0, x0_com, E, order);
    
    % Build Monomial bases for quadrature points (Q) and boundary
    % points (Q0). 
    % -- Draft equation (1)
    % -- Shape matching section 4.3
    Y = monomial_basis_matrix(V', x0_com, order, k);
    Y0 = monomial_basis_matrix(x0, x0_com, order, k);
    
    % For each sampled point, compute the weighting with respect to all
    % shapes.
    % -- Draft Equation (3)
    w = compute_projected_weights(x0, E, V');
    w_x = compute_projected_weights(x0, E, x0);
    [W, ~, W_S] = build_weight_matrix(w, d, k, 'Truncate', false);
    [W0, ~, W0_S] = build_weight_matrix(w_x, d, k, 'Truncate', false);
    
    % Forming gradient of monomial basis with respect to X (undeformed)
    % -- Draft Equation (13)
    dM_dX = monomial_basis_grad(V', x0_com, order);
        
    % Computing each gradient of deformation gradient with respect to
    % projection operator (c are polynomial coefficients)
    dF_dc = vem_dF_dc(dM_dX, W);
    
    % Computing mass matrices
    ME = vem_error_matrix(Y0, W0, W0_S, L);
    M = vem_mass_matrix(Y, W, W_S, L);
    M = sparse((rho*M + k_error*ME));

    ii=1;
    for t=0:dt:30
        b = [];
        for i=1:numel(E)
            b = [b x(:,E{i}) - x0_com];
        end
        b = b(:);
        
        % Solve for polynomial coefficients (projection operators).
        c = L * b;
        
        % Extract deformed center of mass translation.
        p = c(end-d+1:end);
        x_com = x0_com + p;

        % force vector
        dV_dq = zeros(d*(k*numel(E) + 1),1); 
        
        % stiffness matrix
        K = zeros(d*(k*numel(E) + 1), d*(k*numel(E) + 1));
        
        % New deformed positions of quadrature points.
        % (only for visualization)
        Points = zeros(size(V'));
        
        % Computing force and stiffness contributions
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
           
            % Deformation Gradient (Draft equation 13)
            F = dMi_dX * W{i} * W_S{i} * c;
            F = reshape(F,d,d);
            
            % Computing new world position of this point.
            Yi = squeeze(Y(i,:,:))*W{i}*W_S{i}; % weighed monomial basis
            Points(:,i) = Yi * c + x_com;
            
            % Force vector contribution
            dV_dF = neohookean_dF(F,C,D);
            dV_dq = dV_dq + W_S{i}' * dF_dc{i}' * dV_dF;
            
            % Stiffness matrix contribution
            d2V_dF2 = neohookean_dF2(F,C,D);
            K = K - W_S{i}' * dF_dc{i}' * d2V_dF2 * dF_dc{i} * W_S{i};
        end
        K = L' * K * L;
        dV_dq = L' * dV_dq;

        % Error correction force
        x_error = x(:);
        x_error(1:2:end) = x_error(1:2:end) - x_com(1);
        x_error(2:2:end) = x_error(2:2:end) - x_com(2);

        % should just do xx(1:d:end) - x_com?
        f_error = - 2 * ME * x_error;
        f_error = k_error*(dt * P * f_error);
        
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
     
        X_plot.XData = Points(1,:);
        X_plot.YData = Points(2,:);
        drawnow
        
        if save_output
            fn=sprintf('output_png\\100k_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end