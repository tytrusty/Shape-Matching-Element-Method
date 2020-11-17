function vem_sim
    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 17;   	% Lame parameter 1
    D = 0.5 * 150;   	% Lame parameter 2
    gravity = -400;
    k_error = 10000;
    order = 2;
    rho = 100;
    save_output = 0;
    
   % Load mesh
    [V,F] = readOBJ('plane.obj');
    V = V(:,1:2);
    E = boundary_faces(F);
    
    fig=figure(1);
    clf;
    X_plot = plot(V(:,1),V(:,2),'.');
    hold on;
        
    % Get undeformed boundary vertices
    x0 = V(unique([E(:,1) E(:,2)]),:);
%     x0 = x0(:,1:2)';
%     min_I = find(x0(2,:) == max(x0(2,:)));
%     min_I = find(x0(1,:) == max(x0(1,:)));
    min_I = [1 2];
    x0 = [0 0; 0 2;  2 2; 2 0; ]';
    
    % Initial deformed positions and velocities
    x = x0;
    v = zeros(size(x));

    % Create plots.
    E_plot = plot([x(1,:) x(1,1)],[x(2,:) x(2,1)],'o','LineWidth',2, 'Color', 'm');
    axis equal
    % xlim([-1 3])
    ylim([-4 2.5])
    
    % Constraint matrix for boundary conditions.
	% P = fixed_point_constraint_matrix(x0', [2 3]');
    P = fixed_point_constraint_matrix(x0',sort(min_I)');
    
    % Gravity force vector.
  	f_gravity = repmat([0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
    
    % Build Monomial bases
    Q0 = x0 - x0_com;
    if order == 2
        Q_ = zeros(5, size(Q0,2));
        Q_(1:2,:) = Q0;
        Q_(3,:) = Q0(1,:).^2;
        Q_(4,:) = Q0(2,:).^2;
        Q_(5,:) = Q0(1,:).*Q0(2,:);
        Q0=Q_;
    end
    
    % Build 'B' matrices (Aqq in shape matching)
    % Single shape currently.
    [SU, S, SV] = svd(Q0*Q0');
    S = diag(S);
    S = 1./ max(S, 1e-4);
    B = Q0' * (SV * diag(S) * SU');
    
    % Build Monomial bases for all quadrature points
    Q = V' - x0_com;
    if order == 2
        Q_ = zeros(5, size(Q,2));
        Q_(1:2,:) = Q;
        Q_(3,:) = Q(1,:).^2;
        Q_(4,:) = Q(2,:).^2;
        Q_(5,:) = Q(1,:).*Q(2,:);
        Q=Q_;
    end
    
    % Computing inverted mass matrices
    %     M = zeros(numel(x0),numel(x0));
    %     for i = 1:size(Q,2)
    %         BM = B * Q(:,i);
    % %         M = M +  mass_matrix_n(BM);
    %         M = M +  mass_matrix_36(BM);
    %         
    %     end
    %     M = M * rho;
    M = rho * eye(numel(x0));
    
    m = size(Q,2);
    
    % Forming gradient of monomial basis w.r.t X
    if order == 1
       dM_dX = zeros(m,2,2); 
    else
       dM_dX = zeros(m,5,2);  
    end
    for i = 1:m
        factor = (m-1)/m;
        dMi_dX = factor * eye(2);
        if order == 2
            dMi_dX = zeros(5,2); 
            dMi_dX(1:2,:) = factor * eye(2);
            dMi_dX(3,:) = [2*factor*Q(1,i) 0];
            dMi_dX(4,:) = [0 2*factor*Q(2,i)];
            dMi_dX(5,:) = [Q(2,i)*factor Q(1,i)*factor];
        end
        dM_dX(i,:,:) = dMi_dX;
    end
    
    ii=1;
    for t=0:dt:30
        x_com = mean(x,2);
        
        % Compute global shape matching matrix
        p = x - x_com;
        A = p * B;
        %A = A / sqrt(det(A*A')^(1/2)); % preserve volume

        dV_dq = zeros(numel(x),1);     % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        
        % Computing force dV/dq for each point.
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
           
            % Deformation Gradient
            F = A * dMi_dX;
           
            % Force vector
            dV_dF = neohookean_dF(F,C,D);
            dF_dq = vem_dF_dq(B, dMi_dX);
            dV_dq = dV_dq + dF_dq' * dV_dF; % assuming constant area

            % Stiffness matrix
            d2V_dF2 = neohookean_dF2(F,C,D);
            K = K - dF_dq' * d2V_dF2 * dF_dq;
        end
        
        % Error correction force
        f_error = A*Q0 + x_com - x;
        f_error = k_error*(dt * P * f_error(:));
        
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
        
        % Update plot.
        E_plot.XData = [x(1,:) x(1,1)];
        E_plot.YData = [x(2,:) x(2,1)];
        
        % Different 'F' for each x, so this need to be in the loop
        % Points = F *(V'-x0_com) + x_com;
        Points = A * Q + x_com;
        X_plot.XData = Points(1,:);
        X_plot.YData = Points(2,:);
        drawnow
        
        if save_output
            fn=sprintf('output_png\\side_pin_quad_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end