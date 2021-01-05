function vem_nurbs
    % NOTE: gonna deprecate this soon and instead make scripts in examples/
    %       that call vem_simulate_nurbs.m instead.

    % Simulation parameters
    dt = 0.01;          % timestep
%     C = 0.5 * 10000;     % Lame parameter 1
%     D = 0.5 * 150000;    % Lame parameter 2
%     gravity = -100;      % gravity force (direction is -z direction)
%     k_error = 100000;   % stiffness for stability term
%     order = 1;          % (1 or 2) linear or quadratic deformation
%     rho = .5;           % per point density (currently constant)
    C = 0.5 * 5000;     % Lame parameter 1
    D = 0.5 * 150000;    % Lame parameter 2
    gravity = -100;      % gravity force (direction is -z direction)
    k_error = 100000;   % stiffness for stability term
    order = 1;          % (1 or 2) linear or quadratic deformation
    rho = .1;           % per point density (currently constant)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    save_obj = 0;        % (0 or 1) whether to output obj files
    obj_res = 18;       % the amount of subdivision for the output obj

    % Read in NURBs 
    fig=figure(1);
    clf;
    
    % Read in file, generate coordinates and trimesh.
    iges_file = 'rounded_cube.iges';
%     iges_file = 'mug.iges';
%     iges_file = 'cylinder2.iges';
%     iges_file = 'wrench_shrink.igs';
%     iges_file = 'mug3.iges';

    parts=nurbs_from_iges(iges_file);
    parts=nurbs_plot(parts);
    
    % Assembles global generalized coordinates
    [J, q, E, x0] = nurbs_assemble_coords(parts);

    % Initial deformed positions and velocities
    x = x0;
    qdot=zeros(size(q));
    
    % Setup pinned vertices constraint matrix
    [~,I] = mink(x0(1,:),20);
    pin_I = I(1:2);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plot all vertices
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    % Raycasting quadrature as described nowhere yet :)
    [V, vol] = raycast_quadrature(parts, [9 9], 5);

    % Things are hard to tune right now when I produce crazy high volumes
    % so I'm normalizing them (for now...)
    vol = vol ./ max(vol);  

    V=x0;
    vol=ones(size(V,2),1);
    % plot3(V(1,:),V(2,:),V(3,:),'.','Color','r','MarkerSize',20);
    
    % Lame parameters concatenated.
    params = [C, D];
    params = repmat(params,size(V,2),1);
        
    % Gravity force vector.
  	f_gravity = repmat([0 0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
        
    % Shape Matrices
    %E=cell(1);
    %E{1}=1:size(x0,2);
    [B,L] = compute_shape_matrices(x0, x0_com, E, order);
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V, x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order); 
    
    % Compute Shape weights
    a = nurbs_blending_weights(parts, V', 10);
    a_x = nurbs_blending_weights(parts, x0', 10);
    
    % Fixed x values.
    x_fixed = zeros(size(x0));
    for i = 1:numel(pin_I)
        x_fixed(:,pin_I(i))=x0(:,pin_I(i));
    end
    
    % Applying fixed point constraints to NURBS jacobian.
    J = P * J;
    m = size(Q,2);
    
    % Forming gradient of monomial basis w.r.t X
    dM_dX = monomial_basis_grad(V, x0_com, order);
    
    % Computing gradient of deformation gradient w.r.t configuration, q
    % Cover your EYES this code is a disaster
    d = 3;  % dimension (2 or 3)
    dF_dq = vem_dF_dq(B, dM_dX, E, size(x,2), a);
    dF_dq = permute(dF_dq, [2 3 1]);
    dF = cell(m,1);
    dF_I = cell(m,1);
    for i = 1:m
       m1 = dF_dq(:,:,i);
       mm1 = max(abs(m1),[],1);
       I = find(mm1 > 1e-4);
       [~,I] = maxk(mm1,60);
       m1(:,setdiff(1:numel(x),I))=[];
       %sum(mm1 < 1e-4)
       dF_I{i} = I';
       dF{i} = m1;
    end

    % Compute mass matrices
    ME = vem_error_matrix(B, Q0, a_x, d, size(x,2), E);
    M = vem_mass_matrix(B, Q, a, d, size(x,2), E);
    M = ((rho*M + k_error*ME)); %sparse?, doesn't seem to be right now
    % Save & load these matrices for large models to save time.
    %         save('saveM.mat','M');
    %         save('saveME.mat','ME');
    %         M = matfile('saveM.mat').M;
    %         ME = matfile('saveME.mat').ME;

    k=3;
    if order == 2
        k = 9;
    end
    
    xcom_plt = plot3(x0_com(1),x0_com(2),x0_com(3),'.','Color','g','MarkerSize',20);

    ii=1;
    for t=0:dt:30
        tic
        x_com = mean(x,2);
        
        % Compute shape matching matrices
        A=zeros(d, k, numel(E));
        for i=1:numel(E)
            p = x(:,E{i}) - x_com;
            Ai = p*B{i};
            A(:,:,i) = Ai;
        end

        % Preparing input for stiffness matrix mex function.
        n=size(x0,2);
        dF_dqij = permute(dF_dq, [3 1 2]);
%         Aij = permute(A, [3 1 2]);
%         Aij = Aij(:,:);
        
        b = [];
        for i=1:numel(E)
            b = [b x(:,E{i}) - x0_com];
        end
        b = b(:);

        PI = L * b;
        p = PI(end-d+1:end);
        PI = PI(1:end-d);
        PI = reshape(PI, k, d, []);
        PI = permute(PI, [2 1 3]);
        
        norm(PI(:)-A(:))

        Aij = permute(PI, [3 1 2]);
        Aij = Aij(:,:);

        xcom_plt.XData = p(1)+x0_com(1);
        xcom_plt.YData = p(2)+x0_com(2);
        xcom_plt.ZData = p(3)+x0_com(3);
        
        % Stiffness matrix (mex function)
        K = -vem3dmesh_neohookean_dq2(Aij, dF_dqij(:,:), ...
                dM_dX(:,:), a, vol, params,k,n,dF,dF_I);

        % Force vector
        dV_dq = zeros(numel(x),1);

        % Computing force dV/dq for each point.
        % TODO: move this to C++ :)
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
            
            %Aij = zeros(size(A(:,:,1)));
            Aij = zeros(d,k);
            for j = 1:size(E,1)
                %Aij = Aij + A(:,:,j) * a(i,j);
                Aij = Aij + PI(:,:,j) * a(i,j);
            end
            
            % Deformation Gradient
            F = Aij * dMi_dX;
                        
            % Force vector
            dV_dF = neohookean_tet_dF(F,C,D);
            dV_dq = dV_dq + dF_dq(:,:,i)' * dV_dF; % assuming constant area
        end
  
        % Error correction force
        p = x(:);
        p(1:d:end) = p(1:d:end) - x_com(1);
        p(2:d:end) = p(2:d:end) - x_com(2);
        p(3:d:end) = p(3:d:end) - x_com(3);
        f_error = - 2 * ME * p;
        f_error = k_error*(dt * P * f_error(:));
       
        % Force from potential energy.
        f_internal = -dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - dt*dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error);
        qdot = lhs \ rhs;

        % Update position
        q = q + dt*qdot;
        x = reshape(P'*J*q,3,[]) + x_fixed;
        
        % Update NURBs plots
        x_idx=0;
        for i=1:numel(parts)
            x_sz = size(parts{i}.x0,2);
            xi = x(:,x_idx+1:x_idx+x_sz);
            parts{i}.plt.Vertices =xi';
            x_idx = x_idx+x_sz;
        end
        drawnow
        
        if save_obj
            obj_fn = "output/obj/part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,parts,obj_fn,ii);
        end
        
        if save_output
            fn=sprintf('output/img/cylinder_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end
end

