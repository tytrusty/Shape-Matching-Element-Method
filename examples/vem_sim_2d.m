function vem_sim_2d
    % Example 2D VEM simulation.

    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 1700;   	% Lame parameter 1
    D = 0.5 * 15000;   	% Lame parameter 2
    gravity = -100;     % gravity strength
    k_error = 100000;   % Stiffness for VEM stability term
    order = 1;          % (1 or 2) linear or quadratic deformation
    rho = 1;            % density
	d = 2;              % dimension (2 or 3)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    
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
    
    shapes_per_line = 10;
    
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
    [B,L,LM] = compute_shape_matrices(x0, x0_com, E, order);
    
    % Build Monomial bases for quadrature points (Q) and boundary
    % points (Q0). 
    % -- Draft equation (1)
    % -- Shape matching section 4.3
    Q = monomial_basis(V', x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order);     
    
    % For each sampled point, compute the weighting with respect to all
    % shapes.
    % -- Draft Equation (3)
    w = compute_projected_weights(x0, E, V');
    w_x = compute_projected_weights(x0, E, x0);

    % Forming gradient of monomial basis with respect to X (undeformed)
    % -- Draft Equation (13)
    dM_dX = monomial_basis_grad(V', x0_com, order);
        
    % Computing each dF_dq
    dF_dq = vem_dF_dq(B, dM_dX, E, size(x,2), w);
    dF_dq = permute(dF_dq, [2 3 1]);
 
    % Computing mass matrices
    ME = vem_error_matrix(B, Q0, w_x, d, size(x,2), E);
    M = vem_mass_matrix(B, Q, w, d, size(x,2), E);
    M = sparse((rho*M + k_error*ME));

    % The number of elements in the monomial basis.
    % -- example: Draft equation 1 is second order with k == 5
    warning('Need to add truncated weights');
    warning('Move new dm_dx to function');
    k = basis_size(d, order);
    
    % The number of quadrature points where at each point we
    % add stiffness matrix contributions and do all that fun
    % Neohookean stuff.
    m = size(Q,2);
        
    ii=1;
    for t=0:dt:30
        x_com = mean(x,2);
        
        % Compute projection operators for each shape (edge)
        A=zeros(d, k, size(E,1));
        for i=1:size(E,1)
            p = x(:,E{i}) - x_com;
            Ai = p*B{i};
            A(:,:,i) = Ai;
        end
        
        b = [];
        for i=1:numel(E)
            b = [b x(:,E{i}) - x0_com];
        end
        b = b(:);

        PI = L * LM * b;
        p = PI(end-d+1:end);
        PI = PI(1:end-d);
        PI = reshape(PI, k, d, []);
        PI = permute(PI, [2 1 3]);

        dV_dq = zeros(numel(x),1);     % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        K3 = zeros(numel(x), numel(x));
        K2 = zeros(d*(k*numel(E) + 1), d*(k*numel(E) + 1));
        % New deformed positions of quadrature points.
        % (only for visualization)
        Points = zeros(size(V'));
        
        % Computing force and stiffness contributions
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
            
            % Projection operator for this particular point
            % -- Draft equation 3
            Aij = zeros(d,k);
            Wi = zeros(d*k, d*k*numel(E) + d);
            
            
            dM = zeros(d*d, d*k);
            for j = 1:d
                for l = 1:d
                    dM(d*(j-1)+l, k*(l-1)+1:k*l)= dMi_dX(:,j)';
                end
            end
            
            for j = 1:numel(E)
%                 Aij = Aij + A(:,:,j) * w(i,j);
                Aij = Aij + PI(:,:,j) * w(i,j);
                Wi(:, d*k*(j-1)+1:d*k*j) = eye(d*k)*w(i,j);
            end
            
            % Deformation Gradient (Draft equation 13)
            
            F = Aij * dMi_dX;
            
            % Computing new world position of this point.
            Points(:,i) = F * Q(1:d,i) + x0_com + p;
            
            %dfdfqtmp = dF_dq(:,:,i);
            dF_dq_i = dM * Wi; %$ * L * LM;
            %norm(dF_dq_i(:) - dfdfqtmp(:))
            % Force vector contribution
            dV_dF = neohookean_dF(F,C,D);
            dV_dq = dV_dq + dF_dq(:,:,i)' * dV_dF;
            
            % Stiffness matrix contribution
            d2V_dF2 = neohookean_dF2(F,C,D);
            K = K - dF_dq(:,:,i)' * d2V_dF2 * dF_dq(:,:,i);
            K2 = K2 - dF_dq_i' * d2V_dF2 * dF_dq_i;
        end
        LLM = L * LM;
        K = LLM' * K2 * LLM;
        
        % Error correction force
        xx = x(:);
        xx(1:2:end) = xx(1:2:end) - x_com(1);
        xx(2:2:end) = xx(2:2:end) - x_com(2);
        f_error = - 2 * ME * xx;
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