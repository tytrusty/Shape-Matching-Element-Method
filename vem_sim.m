function vem_sim
    dt = 0.01;          % timestep
    C = 0.5 * 1700;            % Lame parameter 1
    D = 0.5 * 15000;            % Lame parameter 2
    
    % M(v-vt) = -dV_dq
    % 1. M formed from shape matching on boundary:
    %   - B fixed for entire shape
    %   - M fixed for each quadrature point
    %   - Need to for P(qdot) in typical manner
    % 2. dV_dq
    %   - dV_dF
    %   - dF_dq
    % 3. Symplectic Euler update
    
   % Load mesh
    [V,F] = readOBJ('plane.obj');
    V = V(:,1:2);
    E = boundary_faces(F);
    
    figure(1)
    clf;
    X_plot = plot(V(:,1),V(:,2),'.');
    hold on;
    
    % Get undeformed boundary vertices
%     x0 = V(unique([E(:,1) E(:,2)]),:);
%     x0 = x0(:,1:2)';
    x0 = [0 0; 0 2;  2 2; 2 0; ]';
    x = x0;
    
    E_plot = plot([x(1,:) x(1,1)],[x(2,:) x(2,1)],'-','LineWidth',2, 'Color', 'm');
    axis equal
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
    plot(x0_com(1), x0_com(2), 'x');
    
    % Build Monomial bases
    order=1;
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
        
    end
    
    M_inv = pinv(rho * M);
    %M_inv2 = inv(rho * M);
    v = zeros(size(x));
    
    for t=0:dt:30
        x_com = mean(x,2);
        
        
        % Compute global shape matching matrix
        P = x - x_com;
        A = P * B;
%         A = A / sqrt(det(A*A')^(1/2)); % preserve volume
        
        % Currently this is just identity, so is constant across the mesh.
        % With higher order we will have a dM_dX for each X. 
        dM_dX = eye(2); % for higher order this should not be identity.
        
        F = A * dM_dX;  % deformation gradient;
        dV_dF = neohookean_dF(F,C,D);
        dF_dq = vem_dF_dq(B, dM_dX);
        dV_dq = 20 * dF_dq' * dV_dF;
                
        a = -M_inv * dV_dq;
        a = reshape(a,2,[]);
        v = v + dt*a + dt*[0 -0.9]';
        
        x = x + dt*v;
        x(:,2:3)= x0(:,2:3);

        
        E_plot.XData = [x(1,:) x(1,1)];
        E_plot.YData = [x(2,:) x(2,1)];
        Points = F * Q + x_com;
        X_plot.XData = Points(1,:);
        X_plot.YData = Points(2,:);
        drawnow
%         x1=x(:)
%         x1=reshape(x1,2,[]);
    end
end