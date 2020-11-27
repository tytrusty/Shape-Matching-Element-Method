function vem_nurbs_old
    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 17000;   	% Lame parameter 1
    D = 0.5 * 1500000;   	% Lame parameter 2
    gravity = -100;
    k_error = 10000;
    order = 2;
    rho = 1;
    save_output = 0;
    
    % Plot NURBs mesh
    function nurbs = plot_srf(nurbs)
            for ii = 1:numel(nurbs)
                xi = reshape(nurbs{ii}.x0, 3, nurbs{ii}.subd(1), nurbs{ii}.subd(2));
                plt = surf(squeeze(xi(1,:,:)),squeeze(xi(2,:,:)),squeeze(xi(3,:,:)),'FaceAlpha',0.8,'EdgeColor','none');
                hold on;
                nurbs{ii}.plt=plt;
            end
    end
    
    % Read in tetmesh
    [V,I] = readNODE([data_dir(), '/meshes_tetgen/Puft/head/coarsest.1.node']);
    [T,~] = readELE([data_dir(), '/meshes_tetgen/Puft/head/coarsest.1.ele']);
    E = boundary_faces(T);
    vol = volume(V,T);
    
    % Swap y&z and flip z axis
    V=V';
    V([2 3],:)=V([3 2],:);
    V(3,:) = -V(3,:);

    % Read in NURBs 
    fig=figure(1);
    clf;
    part=nurbs_from_iges('rounded_cube.iges',6);
    part=plot_srf(part);
    
    % Plot
    xlabel('x'); ylabel('y'); zlabel('z');
    set(gcf,'color','w');
    axis equal
    lighting gouraud;
%     shading interp
    
    lightangle(gca,-15,20)
    zlim([-10 150]);
    x0 = zeros(3,0);
    
    q_size = 0;
    J_size = 0;
    
    % build global position vectors
    for i=1:numel(part)
        idx1=q_size+1;
    	q_size = q_size + size(part{i}.p,1);
        J_size = J_size + size(part{i}.J_flat,2) * size(part{i}.J_flat,3);
        idx2=q_size;
        
        % indices into global configuration vector
        part{i}.idx1=idx1;
        part{i}.idx2=idx2;
    end
    
    q = zeros(q_size,1);
    J = zeros(J_size,q_size);
    J_idx = [0 0];
    for i=1:numel(part)
        q(part{i}.idx1:part{i}.idx2,1) = part{i}.p;
        Ji = part{i}.J_flat(:,:)';
        
        % Block indices
        I1=J_idx(1)+1:J_idx(1)+size(Ji,1);
        I2=J_idx(2)+1:J_idx(2)+size(Ji,2);
        
        % Subsitute block into global NURBs jacobian (J) matrix
        J(I1,I2) = Ji;
        J_idx = J_idx + size(Ji);
    end
    J = sparse(J);
    qdot=zeros(size(q));
    
    for i=1:numel(part)
        idx1=size(x0,2)+1;
        x0 = [x0 part{i}.x0];
        idx2=size(x0,2);
        
        % indices into global configuration vector
        part{i}.idx1=idx1;
        part{i}.idx2=idx2;
    end

    % Initial deformed positions and velocities
    x = x0;
    
    % Setup pinned vertices constraint matrix
    pin_I = find(x0(3,:) > 140);
    pin_I = intersect(find(x0(1,:) < -30), find(x0(2,:) > 55));
    %pin_I = find(x0(3,:) < 65);
    %pin_I = find(x0(1,:) < -40)
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
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
    Q = V - x0_com;
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
    
    M = rho * eye(numel(x0));
    
    % Fixed x values.
    x_fixed = zeros(size(x0));
    for i = 1:numel(pin_I)
        x_fixed(:,pin_I(i))=x0(:,pin_I(i));
    end
    
    J = P * J;
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
             
    % Plot all vertices
%     X_plot=plot3(V(1,:),V(2,:),V(3,:),'.');
%     hold on;
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
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
            dV_dF = SNH_tet_dF(F,C,D);
            dF_dq = vem_dF_dq(B, dMi_dX);
            dV_dq = dV_dq + dF_dq' * dV_dF; % assuming constant area

            % Stiffness matrix
            d2V_dF2 = SNH_tet_dF2(F,C,D);
            K = K - dF_dq' * d2V_dF2 * dF_dq;
        end
        
        % Error correction force
        f_error = A*Q0 + x_com - x;
        f_error = k_error*(dt * P * f_error(:));
        
        % Force from potential energy.
        f_internal = -dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - dt*dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error);
        qdot = lhs \ rhs;   % perform solve

        % Update position
        q = q + dt*qdot;
        x = reshape(P'*J*q,3,[]) + x_fixed;
        
        % Update plot.
%         e_plot.XData = x(1,:);
%         e_plot.YData = x(2,:);
%         e_plot.ZData = x(3,:);
        
        % Different 'F' for each x, so this need to be in the loop
        % Points = F *(V'-x0_com) + x_com;
%         Points = A * Q + x_com;
%         X_plot.XData = Points(1,:);
%         X_plot.YData = Points(2,:);
%         X_plot.ZData = Points(3,:);
        
        % Update NURBs surfaces
        x_idx=0;
        for i=1:numel(part)
            x_sz = size(part{i}.x0,2);
            xi = reshape(x(:,x_idx+1:x_idx+x_sz), 3, part{i}.subd(1), part{i}.subd(2));
            part{i}.plt.XData = squeeze(xi(1,:,:));
            part{i}.plt.YData = squeeze(xi(2,:,:));
            part{i}.plt.ZData = squeeze(xi(3,:,:));
            x_idx = x_idx+x_sz;
        end
        drawnow
        
        if save_output
            fn=sprintf('output_png\\3d_pin_corner_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end

