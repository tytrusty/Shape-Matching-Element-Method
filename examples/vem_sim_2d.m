function vem_sim_2d
    % Example 2D VEM simulation.

    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 170;   	% Lame parameter 1
    D = 0.5 * 1500;   	% Lame parameter 2
    gravity = -40;     % gravity strength
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
    E = boundary_faces(F);
        
    % Get undeformed boundary vertices. 
    V_bnd = unique(E(:));
    x0 = V(V_bnd,:)';

    % Remap edge vertex IDs to IDs into 'x0'
    new_Idx = find(V_bnd);
    [toreplace, bywhat] = ismember(E, V_bnd);
    E(toreplace) = new_Idx(bywhat(toreplace));
    
    % Todo -- only use cells for all edge stuff
    E_cell = cell(size(E,1),1);
    for i =1:size(E,1)
       E_cell{i} = E(i,:); 
    end
    
%     V  = [0 0; 0 2;  2 2; 2 0; 1 1 ];
%     
%     % Simpler cube example.
%     x0  = [0 0; 0 2;  2 2; 2 0; ]';
%     E = [1 2 3; 3 4 1];
%     E_cell = {[1 2 3], [3 4 1]}';
%     E = [1 2; 2 3; 3 4; 4 1];
%     E_cell = {[1 2], [2 3], [3 4], [4 1]}';

%     min_I = find(x0(2,:) == max(x0(2,:))); % pin top side
%     min_I = find(x0(1,:) == min(x0(1,:))); % pin left side
    min_I = [1]; % Pinning the top left corner.

    % Form selection matrices for each shape, allowing us to project
    % out the points that are pinned.
    S = cell(size(E,1),1);
    for i=1:size(E,1)
        S{i} = sparse(zeros(numel(x0), size(E,2)*2));
        for j=1:size(E,2)
            idx = E(i,j);
            S{i}(2*idx-1:2*idx,2*j-1:2*j) = eye(2);
        end
    end
    
    % Initial deformed positions and velocities
    x = [x0 [0 0]'];
    v = zeros(size(x));
        
    % Setup figure
    fig=figure(1);
    clf;
    
    % Plot quadrature points
    X_plot = plot(V(:,1),V(:,2),'.');
    hold on;
    
    % Create plot of boundary edges.
    E_lines = plot([x0(1,E(:,1)); x0(1,E(:,2))], ...
                   [x0(2,E(:,1)); x0(2,E(:,2))],'LineWidth',3);
    hold on;
    
    % Plot pinned vertices
    plot(x0(1,min_I), x0(2,min_I),'.','MarkerSize', 30,'Color','r');
    axis equal
    ylim([-4.5 2.5])
    
    % Constraint matrix for boundary conditions.
    P = fixed_point_constraint_matrix(x',sort(min_I)');
    
    % Gravity force vector.
  	f_gravity = repmat([0 gravity], size(x,2),1)';
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
    [B,~] = compute_shape_matrices(x0, x0_com, E_cell, order);
    
    % Build Monomial bases for quadrature points (Q) and boundary
    % points (Q0). 
    % -- Draft equation (1)
    % -- Shape matching section 4.3
    Q = monomial_basis(V', x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order);     
    
    % For each sampled point, compute the weighting with respect to all
    % shapes.
    % -- Draft Equation (3)
    w = compute_projected_weights(x0, E_cell, V');
    w_x = compute_projected_weights(x0, E_cell, x0);

    % Forming gradient of monomial basis with respect to X (undeformed)
    % -- Draft Equation (13)
    dM_dX = monomial_basis_grad(V', x0_com, order);
        
    % Computing each dF_dq
    dF_dq = vem_dF_dq(B, dM_dX, E_cell, size(x0,2), w);
    dF_dq = permute(dF_dq, [2 3 1]);
 
    % Computing mass matrices
    ME = vem_error_matrix(B, Q0, w_x, d, size(x0,2), E_cell);
    M = vem_mass_matrix(B, Q, w, d, size(x0,2), E_cell);
    M = rho*M + k_error*ME;

    % The number of elements in the monomial basis.
    % -- example: Draft equation 1 is second order with k == 5
    % TODO: compute this the right way :) 
    k=2;
    if order == 2
        k = 5;
    end
    
    % The number of quadrature points where at each point we
    % add stiffness matrix contributions and do all that fun
    % Neohookean stuff.
    m = size(Q,2);
    
    x_com = zeros(d,1);
    com_plt = plot(x0_com(1), x0_com(2), ...
                   '.', 'Color', 'g', 'MarkerSize', 20);
    
    ii=1;
    for t=0:dt:30

        % Compute projection operators for each shape (edge)
        A=zeros(d, k, size(E,1));
        for i=1:size(E,1)
            p = x(:,E(i,:)) - x0_com - x_com;
            Ai = p*B{i};
            A(:,:,i) = Ai;
        end
        
        dV_dq = zeros(numel(x),1);      % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        
        % New deformed positions of quadrature points.
        % (only for visualization)
        Points = zeros(size(V'));
        
        % Computing force and stiffness contributions
        for i = 1:m
            
            % Projection operator for this particular point
            % -- Draft equation 3
            Aij = zeros(size(A(:,:,1)));
            for j = 1:size(E,1)
                Aij = Aij + A(:,:,j) * w(i,j);
            end
            
            % Deformation Gradient (Draft equation 13)
            dMi_dX = squeeze(dM_dX(i,:,:));
            F = Aij * dMi_dX;
            
            % Computing new world position of this point.
            Points(:,i) = F * Q(1:d,i) + x_com + x0_com;
                        
            % Force vector contribution
            dV_dF = neohookean_dF(F,C,D);
            dV_dq = dV_dq + dF_dq(:,:,i)' * dV_dF;

            % Stiffness matrix contribution
            d2V_dF2 = neohookean_dF2(F,C,D);
            K = K - dF_dq(:,:,i)' * d2V_dF2 * dF_dq(:,:,i);
        end
        
        % Error correction force
        xx = x(:);
        xx(1:2:end-3) = xx(1:2:end-3) - x0_com(1);
        xx(2:2:end-2) = xx(2:2:end-2) - x0_com(2);
        f_error = - 2 * ME *  xx;
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
        
        x_com = x(end-1:end)';
        
        % Update plots
        for i = 1:numel(E_lines)
            E_lines(i).XData = x(1,E(i,:)); 
            E_lines(i).YData = x(2,E(i,:)); 
        end
        com_plt.XData = x0_com(1) + x(end-1);
        com_plt.YData = x0_com(2) + x(end);
     
        X_plot.XData = Points(1,:);
        X_plot.YData = Points(2,:);
        drawnow
        
        if save_output
            fn=sprintf('output\\img\\newcom_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end