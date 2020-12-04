function vem_nurbs
    % Simulation parameters
    dt = 0.01;      	% timestep
    C = 0.5 * 17000;   	% Lame parameter 1
    D = 0.5 * 150000;   	% Lame parameter 2
    gravity = -1000;
    k_error = 1000000;
    k_weld = 0;
    order = 1;
    rho = .1;
    save_output = 0;
    save_obj = 0;

    % Read in tetmesh
    [V,I] = readNODE([data_dir(), '/meshes_tetgen/Puft/head/coarsest.1.node']);
    [T,~] = readELE([data_dir(), '/meshes_tetgen/Puft/head/coarsest.1.ele']);
%     [V,I] = readNODE('C:\Users\TY\Desktop\tetgen1.6.0\build\Release\rocket.1.node');
%     [T,~]  = readELE('C:\Users\TY\Desktop\tetgen1.6.0\build\Release\rocket.1.ele');
    E = boundary_faces(T);
    vol = volume(V,T);
       
    % Swap y&z and flip z axis
    V=V';
    V([2 3],:)=V([3 2],:);
    V(3,:) = -V(3,:);

    % Read in NURBs 
    fig=figure(1);
    clf;
        
    %part=nurbs_from_iges('rocket_4.iges',5,0);
%     res=repelem(5,14); res(1)=6;
%     part=nurbs_from_iges('rocket_4.iges',res,0);
    part=nurbs_from_iges('rounded_cube.iges',5,0);
%     res=[8 20];
%     part=nurbs_from_iges('mug.iges',res,1);
    part=plot_nurbs(part);
    
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
        E{i}=idx1:idx2;
    end

    % Initial deformed positions and velocities
    x = x0;
    
    % Setup pinned vertices constraint matrix
    kth_min = mink(x0(3,:),50);
    pin_I = find(x0(3,:) < kth_min(8));
    %     pin_I = find(x0(3,:) < 2);
    %     pin_I = find(x0(1,:) > 6 & x0(3,:) > 6 & x0(3,:) < 8);
    %     pin_I = find(x0(1,:) > max(x0(1,:)) - 1e-4);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plot all vertices
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    %V=[V x0];
%        Vs=x0;
    %plot3(V(1,:),V(2,:),V(3,:),'.');

    % Gravity force vector.
  	f_gravity = repmat([0 0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
        
    % Shape Matrices
    %E=cell(1);
    %E{1}=1:size(x0,2);
    [B,~] = compute_shape_matrices(x0, x0_com, E, order);
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V, x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order); 
    
    % Compute Shape weights
    a = compute_projected_weights(x0, E, V);
    a_x = compute_projected_weights(x0, E, x0);
    
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
    %     S_w = cell(numel(E),1);
    %     S_I = identify_seams(part, x0, E, 3,100);
    %     for i=1:numel(E)
    %         S_w{i} = sparse(zeros(numel(x0), numel(S_I{i})*3));
    %         for j=1:numel(S_I{i})
    %             idx = S_I{i}(j);
    %             S_w{i}(3*idx-2:3*idx,3*j-2:3*j) = eye(3);
    %         end
    %         h1=plot3(x0(1,S_I{i}),x0(2,S_I{i}),x0(3,S_I{i}),'.','Color','m','MarkerSize',30);
    %         delete(h1);
    %     end
    
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
        factor = 1;
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
    
    % Computing gradient of deformation gradient w.r.t configuration, q
    d = 3;  % dimension (2 or 3)
    dF_dq = vem_dF_dq(B, dM_dX, E, size(x,2), a);
    dF_dq = permute(dF_dq, [2 3 1]);
        
%     dFij_dq_sparse = cell(m,1);
%         dFij_dq_sparse{i}=dFij_dq{i};
%         dFij_dq_sparse{i}(dFij_dq_sparse{i} < 1e-8) = 0;
        %colsum=sum(dFij_dq_sparse{i},1);
        %nonz=nnz(colsum);
        %colfind = find(sum(dFij_dq_sparse{i},1) > 0);
        %colsum2=find(sum(dFij_dq_sparse{i},1) > 1e-8);
%         dFij_dq_sparse{i} = sparse(dFij_dq_sparse{i});
%     end

    M = zeros(numel(x), numel(x));
    ME = zeros(numel(x), numel(x));
    for i = 1:size(E,1)
        n=size(B{i},1);
        w_j = reshape(a(:,i)', [1 1 size(Q,2)]);
    	MJ = vem_jacobian(B{i},Q,n,d,size(x,2),E{i});
        MJ = sum(bsxfun(@times, MJ,w_j),3);
        M = M + MJ'*MJ;

        % Stability term
        JE = vem_jacobian(B{i},Q0,n,d,size(x,2),E{i});
        for j=1:size(x0,2)
            I = zeros(d,numel(x0));
            I(:, d*j-2:d*j) = eye(d);
            ME_J = a_x(j,i) * (I - JE(:,:,j));
            ME = ME + ME_J'*ME_J;
        end
    end
    M = ((rho*M + k_error*ME)); %sparse?, doesn't seem to be
%     save('saveM.mat','M');
%     save('saveME.mat','ME');
%     M = matfile('saveM.mat').M;
%     ME = matfile('saveME.mat').ME;
%     M = rho * eye(numel(x0));
    
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
        n=size(x0,2);
        Aij = permute(A, [3 1 2]);
        Aij = Aij(:,:);        
        vol=ones(size(a,1),1);
        params = [C, D];
        params = repmat(params,size(a,1),1);
        
        % Force vector
        % Stiffness matrix
        %K0 = -vem3dmesh_neohookean_dq2(Aij, dFij_dq_orig(:,:), dM_dX(:,:), a, vol, params,k,n,dFij_dq_sparse);

        % Computing force dV/dq for each point.
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
            
            Aij = zeros(size(A(:,:,1)));
            
            for j = 1:size(E,1)
                Aij = Aij + A(:,:,j) * a(i,j);
            end
            
            % Deformation Gradient
            F = Aij * dMi_dX;
                        
            % Force vector
            dV_dF = neohookean_tet_dF(F,C,D);

            dV_dq = dV_dq + dF_dq(:,:,i)' * dV_dF; % assuming constant area

            % Stiffness matrix
            d2V_dF2 = neohookean_tet_dF2(F,C,D);
            K = K - dF_dq(:,:,i)' * d2V_dF2 * dF_dq(:,:,i);
            %K0 = dFij_dq_sparse{i}' * d2V_dF2 * dFij_dq_sparse{i};
        end
        
        % Error correction force
        f_error = - 2 * ME * x(:);
        f_error = k_error*(dt * P * f_error(:));
        
        % Weld force
        f_weld = zeros(numel(x),1);
%         for j = 1:size(E,1)
%             x_pred = A(:,:,j) * Q0(:,S_I{j}) + x_com;
%             x_actual = x(:,S_I{j});
%             diff = x_pred - x_actual;
%             f_weld = f_weld + S_w{j} * diff(:);
%         end
        f_weld = k_weld*(dt * P * f_weld(:));

        % Force from potential energy.
        f_internal = -dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - dt*dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error + f_weld);
        qdot = lhs \ rhs;

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
        
        if save_obj
            obj_fn = "output/obj/part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,part,24,obj_fn);
        end
        
        if save_output
            fn=sprintf('output_png\\img_%03d.png',ii)
            saveas(fig,fn);
        end
        ii=ii+1;
    end
end

