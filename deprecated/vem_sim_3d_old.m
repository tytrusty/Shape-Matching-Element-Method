function vem_sim_3d_old
    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 17;   	% Lame parameter 1
    D = 0.5 * 150;   	% Lame parameter 2
    gravity = -1000;
    k_error = 10000;
    order = 2;
    rho = 10;
    save_output = 0;
    
   % Load mesh
    [V,T,F] = readMESH([data_dir(), '/meshes_mesh/beam.mesh']);
    E = boundary_faces(T);
    
    %triangle areas
    % TODO -- use tetrahedron centroids instead of V.
    vol = volume(V,T);

    fig=figure(1);
    clf;
%     X_plot = plot3(V(:,1),V(:,2),V(:,3),'.');
%     hold on;
        
    % Get undeformed boundary vertices
    V_bnd = unique(E(:));
    x0 = V(V_bnd,:)';
    
    % changem -- remap vertex indices
    new_Idx = find(V_bnd);
    [toreplace, bywhat] = ismember(E, V_bnd);
    E(toreplace) = new_Idx(bywhat(toreplace));
    
%     % Box vertices
%     x0 =[-5 -1 -1; % front bl
%           5 -1 -1; % front br
%           5 -1  1; % front tr
%          -5 -1  1; % front tl
%          -5  1 -1; % back bl
%          -5  1  1; % back tl
%           5  1 -1; % back br
%           5  1  1; % back tr
%     ]';
% 
%     % Box faces
%     E=[ 1 2 3; % Front 1
%         1 3 4; % Front 2
%         1 4 6; % Left 1
%         1 6 5; % Left 2
%         4 3 8; % Top 1
%         8 6 4; % Top 2
%         1 7 2; % Bottom 1
%         1 5 7; % Bottom 2
%         2 7 8; % Right 1
%         2 8 3; % Right 2
%         7 5 6; % Back 1
%         7 6 8; % Back 2
%         ]; 
    
    [x0,E] = readOBJ('beam.obj');
    x0=x0';
    
    min_I = [find(x0(1,:) == min(x0(1,:))) find(x0(1,:) == max(x0(1,:)))] ;
     min_I = find(x0(3,:) == min(x0(3,:)));
    % min_I = [1 4 5 6]; % Pin Left
    % min_I = [1 5];      % Pin left bottom
    % min_I = [1 2];      % Pin front bottom

    % Initial deformed positions and velocities
    x = x0;
    v = zeros(size(x));

    % Create plots.
    E_plot = tsurf(E, x0');
    axis equal
    shading interp
%     lighting gouraud
    lightangle(gca,43,40)
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    % xlim([-1 3])
    ylim([-3 3])
    zlim([-2 1.5])
    xlim([-7 6]);
    view(35,15)
    
    % Constraint matrix for boundary conditions.
    P = fixed_point_constraint_matrix(x0',sort(min_I)');
    
    % Gravity force vector.
  	f_gravity = repmat([0 0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
    
    % Build Monomial bases
    Q0 = x0 - x0_com;
    if order == 2
        Q_ = zeros(9, size(Q0,2));
        Q_(1:3,:) = Q0;
        Q_(4,:) = Q0(1,:).^2;
        Q_(5,:) = Q0(2,:).^2;
        Q_(6,:) = Q0(3,:).^2;
        Q_(7,:) = Q0(1,:).*Q0(2,:);
        Q_(8,:) = Q0(2,:).*Q0(3,:);
        Q_(9,:) = Q0(3,:).*Q0(1,:);
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
        Q_ = zeros(9, size(Q,2));
        Q_(1:3,:) = Q;
        Q_(4,:) = Q(1,:).^2;
        Q_(5,:) = Q(2,:).^2;
        Q_(6,:) = Q(3,:).^2;
        Q_(7,:) = Q(1,:).*Q(2,:);
        Q_(8,:) = Q(2,:).*Q(3,:);
        Q_(9,:) = Q(3,:).*Q(1,:);
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
       dM_dX = zeros(m,3,3); 
    else
       dM_dX = zeros(m,9,3);  
    end
    for i = 1:m
        factor = (m-1)/m;
        dMi_dX = factor * eye(3);
        if order == 2
            dMi_dX = zeros(9,3); 
            dMi_dX(1:3,:) = factor * eye(3);
            dMi_dX(4,:) = [2*factor*Q(1,i) 0 0];
            dMi_dX(5,:) = [0 2*factor*Q(2,i) 0];
            dMi_dX(6,:) = [0 0 2*factor*Q(3,i)];
            dMi_dX(7,:) = [Q(2,i)*factor Q(1,i)*factor 0];
            dMi_dX(8,:) = [0 Q(3,i)*factor Q(2,i)*factor];
            dMi_dX(9,:) = [Q(3,i)*factor 0 Q(1,i)*factor];
        end
        dM_dX(i,:,:) = dMi_dX;
    end
            
    ii=1;
    for t=0:dt:30
        x_com = mean(x,2);
        
        % Compute global shape matching matrix
        p = x - x_com;
        A = p * B;
        %A = A / sqrt(det(A*A')^(1/3)); % preserve volume

        dV_dq = zeros(numel(x),1);     % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        
        % Computing force dV/dq for each point.       
        for i = 1:m
            
            dMi_dX = squeeze(dM_dX(i,:,:));
            % Deformation Gradient
            F = A * dMi_dX;
           
            % Force vector
            dV_dF = neohookean_tet_dF(F,C,D);
            dF_dq = vem_dF_dq(B, dMi_dX);
            dV_dq = dV_dq + dF_dq' * dV_dF; % assuming constant area

            % Stiffness matrix
            d2V_dF2 = neohookean_tet_dF2(F,C,D);
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
        v = reshape(qdot,3,[]);

        % Update position
        x = x + dt*v;
        
        % Update plot.
        E_plot.Vertices = x';
        
        % Different 'F' for each x, so this need to be in the loop
        % Points = F *(V'-x0_com) + x_com;
        Points = A * Q + x_com;
%         X_plot.XData = Points(1,:);
%         X_plot.YData = Points(2,:);
%         X_plot.ZData = Points(3,:);
        drawnow
        
        if save_output
            fn=sprintf('output_png\\3d_pin_sides_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end