function vem_sim
    dt = 0.01;      	% timestep
    C = 0.5 * 17;   	% Lame parameter 1
    D = 0.5 * 150;   	% Lame parameter 2
    gravity = -100;
    order = 2;          
    
   % Load mesh
    [V,F] = readOBJ('plane.obj');
    V = V(:,1:2);
    E = boundary_faces(F);
    
    fig=figure(1)
    clf;
    X_plot = plot(V(:,1),V(:,2),'.');
    hold on;
        
    % Get undeformed boundary vertices
    x0 = V(unique([E(:,1) E(:,2)]),:);
    x0 = x0(:,1:2)';
    min_I = find(x0(2,:) == max(x0(2,:)));
    x0 = [0 0; 0 2;  2 2; 2 0; ]';
    x = x0;
    
    g_force = repmat([0 gravity], size(x0,2),1)';
    g_force = g_force(:);
    
    E_plot = plot([x(1,:) x(1,1)],[x(2,:) x(2,1)],'o','LineWidth',2, 'Color', 'm');
    axis equal
    xlim([-1 3])
    ylim([-.5 2.5])
    
%    P = fixed_point_constraint_matrix(x0', [2 3]')
%      b = qt - P'*P*qt;
%     qt = P*qt;
%     vt = 0*qt;
%     
    % Undeformed Center of mass
    x0_com = mean(x0,2);
    plot(x0_com(1), x0_com(2), 'x');
    
    % Build Monomial bases
    Q = x0 - x0_com;
    if order == 2
        Q_ = zeros(5, size(Q,2));
        Q_(1:2,:) = Q;
        Q_(3,:) = Q(1,:).^2;
        Q_(4,:) = Q(2,:).^2;
        Q_(5,:) = Q(1,:).*Q(2,:);
        Q=Q_;
    end
    
    % Build 'B' matrices (Aqq in shape matching)
    % Single shape currently.
    [SU, S, SV] = svd(Q*Q');
    S = diag(S);
    S = 1./ max(S, 1e-4);
    B = Q' * (SV * diag(S) * SU');
    

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
    M = zeros(numel(x0),numel(x0));
    rho = 100;
    for i = 1:size(Q,2)
        BM = B * Q(:,i);
        M = M +  mass_matrix_n(BM);
%         M = M +  mass_matrix_36(BM);
        
    end
    M = M * rho;
    
    M_inv = pinv(rho * M);
    %M_inv2 = inv(rho * M);
    v = zeros(size(x));
    ii=1;
    
    M = 1e3 * eye(8);
%     M_inv = 1e-3 * eye(8);
    %project mass matrix down
    %M_inv = P*M_inv*P';
    
    for t=0:dt:30
        x_com = mean(x,2);
        
        % Compute global shape matching matrix
        p = x - x_com;
        A = p * B;
        % A = A / sqrt(det(A*A')^(1/2)); % preserve volume

        dV_dq = zeros(numel(x),1);     % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        
        % Computing force dV/dq for each point.
        m = size(Q,2);
        for i = 1:m
            % Forming gradient of monomial basis w.r.t X
        	factor = (m-1)/m;
            dM_dX = factor * eye(2);
            if order == 2
                dM_dX = zeros(5,2);
                dM_dX(1:2,:) = factor * eye(2);
                dM_dX(3,:) = [2*factor*V(i,1) 0];
                dM_dX(4,:) = [0 2*factor*V(i,2)];
                dM_dX(5,:) = [V(i,2)*factor V(i,1)*factor];
            end
           
            % Deformation Gradient
            F = A * dM_dX;
           
            % Force vector
            dV_dF = neohookean_dF(F,C,D);
            dF_dq = vem_dF_dq(B, dM_dX);
            dV_dq = dV_dq + dF_dq' * dV_dF; % assuming constant area

            % Stiffness matrix
            d2V_dF2 = neohookean_dF2(F,C,D);
            K = K - dF_dq' * d2V_dF2 * dF_dq;
                      
        end
               
%         a = -M_inv * dV_dq;
%         a = reshape(a,2,[]);
%         v = v + dt*a + dt*[0 gravity]';
%         x = x + dt*v;
        
        lhs = (M - dt*dt*K);
        rhs = M*v(:) - dt*dV_dq + dt*g_force;
        qdot = lhs \ rhs;
        v = reshape(qdot,2,[]);
        x = x + dt*v;
        
        x(:,2:3)= x0(:,2:3);
        % x(:,min_I)= x0(:,min_I);
        
        E_plot.XData = [x(1,:) x(1,1)];
        E_plot.YData = [x(2,:) x(2,1)];
        % Different 'F' for each x, so this need to be in the loop
        % Points = F *(V'-x0_com) + x_com;
        Points = A * Q + x_com;
        X_plot.XData = Points(1,:);
        X_plot.YData = Points(2,:);
        drawnow
        
%         fn=sprintf('output_png\\lin_%03d.png',ii)
%         saveas(fig,fn);
%         ii=ii+1;
    end
end