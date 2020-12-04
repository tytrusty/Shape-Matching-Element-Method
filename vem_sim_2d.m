function vem_sim_2d
    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 1700;   	% Lame parameter 1
    D = 0.5 * 150000;   	% Lame parameter 2
    gravity = -1000;
    k_error = 1000000;
    order = 1;
    rho = 10;
    save_output = 0;
    
   % Load mesh
    [V,F] = readOBJ('plane.obj');
    V = V(:,1:2);
    E = boundary_faces(F);
    
    fig=figure(1);
    clf;
    
%     V=[0 0; 0 2; 2 0; 2 2; 1 1];
    X_plot = plot(V(:,1),V(:,2),'.');
    hold on;
        
    % Get undeformed boundary vertices
    V_bnd = unique(E(:));
    x0 = V(V_bnd,:)';
    
    % changem -- remap vertex indices
    new_Idx = find(V_bnd);
    [toreplace, bywhat] = ismember(E, V_bnd);
    E(toreplace) = new_Idx(bywhat(toreplace));
    E_cell = cell(size(E,1),1);
    for i =1:size(E,1)
       E_cell{i} = E(i,:); 
    end
    min_I = find(x0(2,:) == max(x0(2,:)));
    min_I = find(x0(1,:) == min(x0(1,:)));
    min_I = [1];
    %%%% NEW %%%%
%     min_I = [3 4];
%     x0 = [0 0;0 2; 0 2; 2 2; 2 2; 2 0; 2 0; 0 0]';
%     E = [1 2; 3 4; 5 6; 7 8;];
%     E_cell = {[1 2], [3 4], [5 6], [7 8]};

%     E = [7 1 2; 1 3 4; 3 5 6; 5 7 8;];
%     E_cell = {[7 1 2], [1 3 4], [3 5 6], [5 7 8]};    
    %%%% ORIG %%%%
% 	min_I = [2 3];
%     x0 = [0 0; 0 2;  2 2; 2 0; ]';
%     E = [2 1; 3 2; 4 3; 1 4];
%     E_cell = {[1 2], [2 3], [3 4], [4 1]};

    % Form selection matrices for each shape.
    S = cell(size(E,1),1);
    for i=1:size(E,1)
        S{i} = sparse(zeros(numel(x0), size(E,2)*2));
        for j=1:size(E,2)
            idx = E(i,j);
            S{i}(2*idx-1:2*idx,2*j-1:2*j) = eye(2);
        end
    end
    
    % Midpoints of each shape
    E_centroids = zeros(size(V,2), size(E,1));
    for i = 1:size(E,1)
        E_centroids(:,i) = mean(x0(:,E(i,:)),2);
    end
    
    % Initial deformed positions and velocities
    x = x0;
    v = zeros(size(x));

    % Create plots.
%     E_plot = plot([x(1,:) x(1,1)],[x(2,:) x(2,1)],'o','LineWidth',3, 'Color', 'm');
    E_lines = plot([x0(1,E(:,1)); x0(1,E(:,2))], [x0(2,E(:,1)); x0(2,E(:,2))],'LineWidth',3);

    axis equal
%     xlim([-3 3])
    ylim([-4.5 2.5])
    hold on;
    plot(x0(1,min_I), x0(2,min_I),'.','MarkerSize', 30,'Color','r');
    
    % Constraint matrix for boundary conditions.
    P = fixed_point_constraint_matrix(x0',sort(min_I)');
    
    % Gravity force vector.
  	f_gravity = repmat([0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
    
    % Build Monomial bases
    [B,~] = compute_shape_matrices(x0, x0_com, E_cell, order);
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V', x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order);     
    
    m = size(Q,2);
    
    % Compute per-element shape weights.
    %     a = compute_weights(E_centroids,E,V');
    a = compute_projected_weights(x0, E_cell, V');
    a_x = compute_projected_weights(x0, E_cell, x0);
    
    % Forming gradient of monomial basis w.r.t X
    dM_dX = monomial_basis_grad(V', x0_com, order);
        
    % Computing each dF_dq
    d = 2;              % dimension (2 or 3)
    dF_dq = vem_dF_dq(B, dM_dX, E_cell, size(x,2), a);
    dF_dq = permute(dF_dq, [2 3 1]);
 
    M = zeros(numel(x), numel(x));
    ME = zeros(numel(x), numel(x));
    for i = 1:size(E,1)
        n=size(B{i},1);
        w_j = reshape(a(:,i)', [1 1 size(Q,2)]);
    	MJ = vem_jacobian(B{i},Q,n,d,size(x,2),E_cell{i});
        MJ = sum(bsxfun(@times, MJ,w_j),3);
        M = M + MJ'*MJ;
        
        % Stability term
        JE = vem_jacobian(B{i},Q0,n,d,size(x,2),E_cell{i});
        for j=1:size(x0,2)
            I = zeros(d,numel(x0));
            I(:, d*j-1:d*j) = eye(d);
            ME_J = a_x(j,i) * (I - JE(:,:,j));
            ME = ME + ME_J'*ME_J;
        end
    end
    M = sparse((rho*M + k_error*ME));
    %     M = 1000*eye(numel(x));
    k=2;
    if order == 2
        k = 5;
    end
    
    ii=1;
    for t=0:dt:30
        x_com = mean(x,2);
        
        % Compute global shape matching matrices
        A=zeros(d, k, size(E,1));
        for i=1:size(E,1)
            p = x(:,E(i,:)) - x_com;
            Ai = p*B{i};
            A(:,:,i) = Ai;
        end
        
        dV_dq = zeros(numel(x),1);     % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        
        Points = zeros(size(V'));
        
        % Computing force dV/dq for each point.
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
            
            Aij = zeros(size(A(:,:,1)));
            dFij_dq = zeros(4,numel(x));
            for j = 1:size(E,1)
                Aij = Aij + A(:,:,j) * a(i,j);
%                 dFij_dq = dFij_dq + dF_dq(:,:,i,j) * a(i,j) * S{j}';
            end
            
            % Deformation Gradient
            F = Aij * dMi_dX;

            Points(:,i) = F * Q(1:2,i) + x_com;
                        
            % Force vector
            dV_dF = neohookean_dF(F,C,D);

            dV_dq = dV_dq + dF_dq(:,:,i)' * dV_dF; % assuming constant area

            % Stiffness matrix
            d2V_dF2 = neohookean_dF2(F,C,D);
            K = K - dF_dq(:,:,i)' * d2V_dF2 * dF_dq(:,:,i);
        end
        
        % Error correction force
        f_error = - 2 * ME * x(:);
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
            E_lines(i).XData = x(1,E(i,:)); 
            E_lines(i).YData = x(2,E(i,:)); 
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