function test_nurbs_cube_animation
    % Simulation parameters
    dt = 0.005;      	% timestep
    C = 0.5 * 17000;   	% Lame parameter 1
    D = 0.5 * 150000;    % Lame parameter 2
    gravity = -450;
    k_error = 100000;
    order = 2;
    rho = 10;
    save_output = 1;
    save_obj = 0;
    obj_res = 18;
    
    plot_mode = 2;
    plot_modes = [repelem(0,100) repelem(1,600) repelem(2,500)];
    patch_step = [repelem(0,100) repelem(1,100) repelem(2,100) repelem(3,100) repelem(4,100) repelem(5,100) repelem(6,100)];
    % TODO --- use volume of quadrature points ... vol = volume(V,T);

    % Read in NURBs 
    fig=figure(1);
    clf;
    
    
    % Some files I test on
    iges_file = 'rounded_cube.iges';

    % Resolution indicates how many point samples we will take on each
    % e.g. 6 means we have 6 samples in both the U & V coordinates, so
    %      a total of 36 samples across the NURBs patch.
    resolution = 6;
    part=nurbs_from_iges(iges_file, resolution,0);
    part=nurbs_plot(part);
    
    % Raycasting quadrature as described nowhere yet :)
    V = raycast_quadrature(part, [6 6], 5)';
    % plot3(V(1,:),V(2,:),V(3,:),'.','Color','r','MarkerSize',20);

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
    [~,I] = mink(x0(3,:),20);
    pin_I = I(1:4);
    % pin_I = find(x0(1,:) < -2.3 & x0(3,:) > 14 );
    % pin_I = find(x0(1,:) > 6 & x0(3,:) > 6 & x0(3,:) < 8);
    % pin_I = find(x0(1,:) > max(x0(1,:)) - 1e-4);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plot all vertices
%     X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
%     hold on;
    % V=x0;

    % Gravity force vector.
  	f_gravity = repmat([0 0 gravity], size(x0,2),1)';
    f_gravity = dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
         
    [V_cube,F_cube] = readOBJ('models/cube.obj');
    V_cube = V_cube / 2.5;
    V_cubes = cell(numel(part),1);
    Q_cubes = cell(numel(part),1);
    grid off;
    view(-16,7);
    zlim([60 145]);
    set(gca,'visible','off');
    
    part{1}.plt.FaceColor = [1 0 1];
    part{2}.plt.FaceColor = [0 0 1];
    part{3}.plt.FaceColor = [1 0 0];
    part{4}.plt.FaceColor = [0 1 0];
    part{5}.plt.FaceColor = [1 1 0];
    part{6}.plt.FaceColor = [1 0.5 0];
    patch_i = 1;
    for i = 1:numel(part)
        t1 = 0.2;
        center = t1 * mean(part{i}.x0, 2) + (1-t1)*x0_com ;
        dims = max(part{i}.x0,[],2) - min(part{i}.x0,[],2);
        dims = sort(dims);
        median = dims(2);
        V_cubes{i} = V_cube *  0.9;
        V_cubes{i} = V_cubes{i} * median + center';
        Q_cubes{i} = monomial_basis(V_cubes{i}', x0_com, order);   
        part{i}.plt.DiffuseStrength = 1;
        part{i}.plt.EdgeColor = 'none';
    end
    cube_plt = trisurf(F_cube, V_cubes{patch_i}(:,1),V_cubes{patch_i}(:,2),V_cubes{patch_i}(:,3), ...
                      'FaceColor', part{patch_i}.plt.FaceColor, ...
                      'FaceAlpha',0.8,'EdgeColor','none');
    cube_plt.AmbientStrength = 0.5;
    
    for i = [1:patch_i-1 patch_i+1:numel(part)]
        part{i}.plt.FaceAlpha = 0.2;
    end
        
    %plot3(V(1,44),V(2,44),V(3,44),'.','MarkerSize',20);
    
    
    % I don't know how the hell color blending works so I'm using
    % https://colordesigner.io/color-mixer
    % for this.
    
    % Shape Matrices
    %E=cell(1);
    %E{1}=1:size(x0,2);
    [B,~] = compute_shape_matrices(x0, x0_com, E, order);
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V, x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order); 
    
    % Compute Shape weights
    a = nurbs_blending_weights(part, V', 40);
    a_x = nurbs_blending_weights(part, x0', 40);
    
    pnt_w = a(44,:);
    color = (pnt_w(1)*part{1}.plt.FaceColor + ...
            pnt_w(2)*part{2}.plt.FaceColor + ...
            pnt_w(3)*part{3}.plt.FaceColor + ...
            pnt_w(4)*part{4}.plt.FaceColor + ...
            pnt_w(5)*part{5}.plt.FaceColor + ...
            pnt_w(6)*part{6}.plt.FaceColor);
    center = V(:,44);
    V_cube_pnt = V_cube *  0.9;
    V_cube_pnt = V_cube_pnt * 30 + center';
    Q_cube_pnt = monomial_basis(V_cube_pnt', x0_com, order);    
        
    % Form selection matrices for each shape.
    S = cell(numel(E),1);
    for i=1:size(E,1)
        S{i} = sparse(zeros(numel(x0), numel(E{i})*3));
        for j=1:numel(E{i})
            idx = E{i}(j);
            S{i}(3*idx-2:3*idx,3*j-2:3*j) = eye(3);
        end
    end
    
    % Fixed x values.
    x_fixed = zeros(size(x0));
    for i = 1:numel(pin_I)
        x_fixed(:,pin_I(i))=x0(:,pin_I(i));
    end
    
    J = P * J;
    m = size(Q,2);
    
    % Forming gradient of monomial basis w.r.t X
    dM_dX = monomial_basis_grad(V, x0_com, order);
    
    % Computing gradient of deformation gradient w.r.t configuration, q
    % Cover your EYES this code is a disaster
    d = 3;  % dimension (2 or 3)
    dF_dq = vem_dF_dq(B, dM_dX, E, size(x,2), a);
    dF_dq = permute(dF_dq, [2 3 1]);
    SdF = cell(m,1);
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
       SdF{i} = sparse(zeros(numel(I), numel(x)));
       ind=sub2ind(size(SdF{i}), 1:numel(I),I);
       SdF{i}(ind)=1;
       dF{i} = m1;
    end

    % Compute mass matrices
    ME = vem_error_matrix(B, Q0, a_x, d, size(x,2), E);
    M = vem_mass_matrix(B, Q, a, d, size(x,2), E);
    M = ((rho*M + k_error*ME)); %sparse?, doesn't seem to be right now
    %     save('saveM.mat','M');
    %     save('saveME.mat','ME');
    % 	M = matfile('saveM.mat').M;
    % 	ME = matfile('saveME.mat').ME;

    k=3;
    if order == 2
        k = 9;
    end
        
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
        
        plot_mode = plot_modes(ii);
        if plot_mode == 0
            for i = [1:patch_i-1 patch_i+1:numel(part)]
                part{i}.plt.FaceAlpha = 0.9;
            end
            cube_plt.FaceAlpha = 0;
        elseif (plot_mode == 1)
            patch_i = patch_step(ii);
            for i = [1:patch_i-1 patch_i+1:numel(part)]
                part{i}.plt.FaceAlpha = 0.05;
            end
            part{patch_i}.plt.FaceAlpha = 0.65;
            cube_plt.FaceAlpha = 0.8;
            A_cube = A(:,:,patch_i);
            V_cube = A_cube * Q_cubes{patch_i} + x_com;
            cube_plt.Vertices = V_cube';
            cube_plt.FaceColor = part{patch_i}.plt.FaceColor;
        else
            for i = 1:numel(part)
                part{i}.plt.FaceAlpha = 0.1;
            end
            cube_plt.FaceAlpha = 0.95;
            A_cube=zeros(d, k);
            for i=1:numel(E)
                A_cube = A_cube + pnt_w(i)*A(:,:,i);
            end
            V_cube = A_cube * Q_cube_pnt + x_com;
            cube_plt.Vertices = V_cube';
            cube_plt.FaceColor = color;
        end
%         A_cube = A(:,:,patch_i);
%         A_cube=zeros(d, k);
%         for i=1:numel(E)
%             A_cube = A_cube + pnt_w(i)*A(:,:,i);
%         end
%         V_cube = A_cube * Q_cube_pnt + x_com;
%         cube_plt.Vertices = V_cube';
%         cube_plt.FaceColor = color;
                
        % Computing force dV/dq for each point.
        n=size(x0,2);
        dF_dqij = permute(dF_dq, [3 1 2]);
        Aij = permute(A, [3 1 2]);
        Aij = Aij(:,:);        
        vol=ones(size(a,1),1);
        params = [C, D];
        params = repmat(params,size(a,1),1);
        
        % Stiffness matrix
        K = -vem3dmesh_neohookean_dq2(Aij, dF_dqij(:,:), dM_dX(:,:), a, vol, params,k,n,dF,dF_I);
        
        % Force vector
        dV_dq = zeros(numel(x),1);

        % Computing force dV/dq for each point.
        % TODO -- move this to C++ :)
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
        end
  
        % Error correction force
        f_error = - 2 * ME * x(:);
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
            nurbs_write_obj(q,part,obj_res,obj_fn,ii);
        end
        
        if save_output
            fn=sprintf('output/img/fix_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end
end

