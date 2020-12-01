function vem_nurbs
    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 17000;   	% Lame parameter 1
    D = 0.5 * 150000;   	% Lame parameter 2
    gravity = -100;
    k_error = 10000;
    k_weld = 10000;
    order = 2;
    rho = 1;
    save_output = 0;
    
    % Plot NURBs mesh
    function nurbs = plot_srf(nurbs)
            cm=jet(numel(nurbs));
            for ii = 1:numel(nurbs)
                
                xi = reshape(nurbs{ii}.x0, 3, nurbs{ii}.subd(1), nurbs{ii}.subd(2));
                plt = surf(squeeze(xi(1,:,:)),squeeze(xi(2,:,:)),squeeze(xi(3,:,:)),'FaceAlpha',.9,'EdgeColor','black','FaceColor',cm(ii,:));
                hold on;
                nurbs{ii}.plt=plt;
            end
    end
    
    % Read in tetmesh
    [V,I] = readNODE([data_dir(), '/meshes_tetgen/Puft/head/coarsest.1.node']);
    [T,~] = readELE([data_dir(), '/meshes_tetgen/Puft/head/coarsest.1.ele']);
%     [V,I] = readNODE('C:\Users\TY\Desktop\vem\tetgen1.6.0\build\Release\rocket.1.node');
%     [T,~]  = readELE('C:\Users\TY\Desktop\vem\tetgen1.6.0\build\Release\rocket.1.ele');
    E = boundary_faces(T);
    vol = volume(V,T);
       
    % Swap y&z and flip z axis
    V=V';
    V([2 3],:)=V([3 2],:);
    V(3,:) = -V(3,:);
%     V = V/20;

    % Read in NURBs 
    fig=figure(1);
    clf;
        
%     part=nurbs_from_iges('rocket_4.iges',7,0);
%     res=repelem(7,14);
%     res(5)=7; res(1)=11;
%     part=nurbs_from_iges('rocket_4.iges',res,0);
%     part=nurbs_from_iges('rocket.iges',res,0);
    part=nurbs_from_iges('cylinder.iges',14,0);
%     part=nurbs_from_iges('rounded_cube.iges',5,0);
    part=plot_srf(part);
    
    % Plot

    set(gcf,'color','w');
    axis equal
    lighting gouraud;
%     shading interp
    
    lightangle(gca,72,-4)
%     zlim([0 90]);
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
    E=cell(numel(part),1);
    
    for i=1:numel(part)
        idx1=size(x0,2)+1;
        x0 = [x0 part{i}.x0];
        idx2=size(x0,2);
        
        % indices into global configuration vector
        part{i}.idx1=idx1;
        part{i}.idx2=idx2;
        E{i}=idx1:idx2;
    end
    
    %V=[V x0];
    V=x0;
%     plot3(V(1,:),V(2,:),V(3,:),'.');

    % Initial deformed positions and velocities
    x = x0;
    
    % Setup pinned vertices constraint matrix
    kth_min = mink(x0(3,:),8);
    pin_I = find(x0(3,:) < kth_min(8));
    %pin_I = intersect(find(x0(1,:) < -30), find(x0(2,:) > 55));
    %pin_I = find(x0(3,:) < 5);
    %pin_I = find(x0(1,:) < -40)
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Gravity force vector.
  	f_gravity = repmat([0 0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
        
    % Shape Matrices
    %E=cell(1);
    %E{1}=1:size(x0,2);
    [B,Q0] = compute_shape_matrices(x0, x0_com, E, order);
    Q0_stacked = [];
    for i = 1:numel(Q0)
        Q0_stacked = [Q0_stacked Q0{i}];
    end
    
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V, x0_com, order);
    
    % Compute Shape weights
    a = compute_projected_weights(x0, E, V);
    
    % Form selection matrices for each shape.
    S = cell(numel(E),1);
    for i=1:size(E,1)
        S{i} = sparse(zeros(numel(x0), numel(E{i})*3));
        for j=1:numel(E{i})
            idx = E{i}(j);
            S{i}(3*idx-2:3*idx,3*j-2:3*j) = eye(3);
        end
    end
    
    % Form selection matrices for each shape.
    S_w = cell(numel(E),1);
    S_I = identify_seams(part, x0, E,1,100);
    for i=1:numel(E)
        S_w{i} = sparse(zeros(numel(x0), numel(S_I{i})*3));
        for j=1:numel(S_I{i})
            idx = S_I{i}(j);
            S_w{i}(3*idx-2:3*idx,3*j-2:3*j) = eye(3);
        end
        h1=plot3(x0(1,S_I{i}),x0(2,S_I{i}),x0(3,S_I{i}),'.','Color','m','MarkerSize',30);
        delete(h1);
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
        factor = 1;%(m-1)/m;
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
    
    % Computing each dF_dq
    d = 3;  % dimension (2 or 3)
    dF_dq = cell(numel(E),1);
    % Init cell matrices
    for i = 1:size(E,1)
        n=size(B{i},1);
        dF_dq{i} = zeros(d*d,d*n,m);
    end
    for i = 1:m
        dMi_dX = squeeze(dM_dX(i,:,:));
        % Per-shape forces & stiffness contribution.
        for j = 1:size(E,1)
            dF_dq{j}(:,:,i) = vem_dF_dq(B{j}, dMi_dX);
        end
    end
    
    
    % Plot all vertices
%     X_plot=plot3(V(1,:),V(2,:),V(3,:),'.');
%     hold on;
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    k=3;
    if order == 2
        k = 9;
    end
        
    ii=1;
    for t=0:dt:30
        x_com = mean(x,2);
        
        % Compute shape matching matrices
        A=zeros(d, k, numel(E));
        for i=1:numel(E)
            p = x(:,E{i}) - x_com;
            Ai = p*B{i};
            A(:,:,i) = Ai;
        end
        
        dV_dq = zeros(numel(x),1);     % force vector
        K = zeros(numel(x), numel(x)); % stiffness matrix
        
        % Computing force dV/dq for each point.
        dM_dX_flat = dM_dX(:,:);
        for j = 1:size(E,1)
            Aj = A(:,:,j);
            dFj_dq = permute(dF_dq{j}, [3 1 2]);
            dFj_dq=dFj_dq(:,:);
            w = a(:,j);
            vol=ones(size(w));
            params = [C, D];
            params = repmat(params,numel(w),1);
            n=size(B{j},1);
            
            % Force vector
            dV_dq = dV_dq + S{j} * vem3dmesh_neohookean_dq(Aj, dFj_dq, dM_dX_flat, w, vol, params,k,n);
            % Stiffness matrix
            d2V_q2 = vem3dmesh_neohookean_dq2(Aj, dFj_dq, dM_dX_flat, w, vol, params,k,n);
            K = K - S{j} * d2V_q2 * S{j}';
        end
        
        % Error correction force
        f_error = zeros(numel(x),1);
        for j = 1:size(E,1)
            x_pred = A(:,:,j) * Q0{j} + x_com;
            x_actual = x(:,E{j});
            diff = x_pred - x_actual;
            f_error = f_error + S{j} * diff(:);
        end
        f_error = k_error*(dt * P * f_error(:));
        
        % Weld force
        f_weld = zeros(numel(x),1);
        for j = 1:size(E,1)
            x_pred = A(:,:,j) * Q0_stacked(:,S_I{j}) + x_com;
            x_actual = x(:,S_I{j});
            diff = x_pred - x_actual;
            f_weld = f_weld + S_w{j} * diff(:);
        end
        f_weld = k_weld*(dt * P * f_weld(:));
        
        
        % Force from potential energy.
        f_internal = -dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - dt*dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error + f_weld);
        qdot = lhs \ rhs;   % perform solve

        % Update position
        q = q + dt*qdot;
        x = reshape(P'*J*q,3,[]) + x_fixed;
        
        % Update NURBs plots
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
            fn=sprintf('output_png\\pin_bottom_fail_quad_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end

