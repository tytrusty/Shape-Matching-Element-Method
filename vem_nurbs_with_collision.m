function vem_nurbs
    % Simulation parameters
    dt = 0.01;          % timestep
    C = 0.5 * 1700;    % Lame parameter 1
    D = 0.5 * 15000;   % Lame parameter 2
    gravity = -50;     % gravity force (direction is -z direction)
    k_error = 10000;   % stiffness for stability term
    order = 1;          % (1 or 2) linear or quadratic deformation
    rho = .1;            % per point density (currently constant)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    save_obj = 0;       % (0 or 1) whether to output obj files
    obj_res = 18;       % the amount of subdivision for the output obj

    % Read in NURBs 
    fig=figure(1);
    clf;
    
    % Some files I test on
    iges_file = 'rounded_cube.iges';
    % iges_file = 'starship_nose.iges';
    % iges_file = 'puft_simple.iges';
    % iges_file = 'castle_simple.iges';
    % iges_file = 'rocket_with_nose.iges';
    % iges_file = 'mug.iges';
    % iges_file = 'rocket.iges'; % this rocket doesn't have the nosecone
    
    % Resolution indicates how many point samples we will take on each
    % e.g. 6 means we have 6 samples in both the U & V coordinates, so
    %      a total of 36 samples across the NURBs patch.
    resolution = 6;
    
    % NOTE: if you a warning saying "Matrix is singular", this means the
    %       resolution is too low. I will fix this so that resolution is
    %       set automatically on monday :)

    % resolution = repelem(8,9); resolution(1)=11; % resolution(17)=11;
    part=nurbs_from_iges(iges_file, resolution, 0);
    part=nurbs_plot(part);
    
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
    [~,I] = mink(x0(1,:),20);
    pin_I = I(1:9);
    % pin_I = find(x0(1,:) < -2.3 & x0(3,:) > 14 );
    % pin_I = find(x0(1,:) > 6 & x0(3,:) > 6 & x0(3,:) < 7);
    % pin_I = find(x0(1,:) == max(x0(1,:)));
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plot all vertices
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    % Raycasting quadrature as described nowhere yet :)
    [V, vol] = raycast_quadrature(part, [9 9], 5);

    % Things are hard to tune right now when I produce crazy high volumes
    % so I'm normalizing them (for now...)
    vol = vol ./ max(vol);  

    % plot3(V(1,:),V(2,:),V(3,:),'.','Color','r','MarkerSize',20);
    %V=x0;
    %vol=ones(size(V,2),1);
    
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
    [B,~] = compute_shape_matrices(x0, x0_com, E, order);
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V, x0_com, order);
    Q0 = monomial_basis(x0, x0_com, order); 
    
    % Compute Shape weights
    a = nurbs_blending_weights(part, V', 30);
    a_x = nurbs_blending_weights(part, x0', 30);
    
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
    % Save & load these matrices for large models to save time.
    % save('saveM.mat','M');
    % save('saveME.mat','ME');
    % M = matfile('saveM.mat').M;
    % ME = matfile('saveME.mat').ME;

    k=4;
    if order == 2
        k = 10;
    end
    
    % Sphere for collision
    sphere_c = [0;0;15];
    sphere_r = 1;  
    [sphere_x, sphere_y, sphere_z] = sphere;
    sphere_x = sphere_x * sphere_r + sphere_c(1);
    sphere_y = sphere_y * sphere_r + sphere_c(2);
    sphere_z = sphere_z * sphere_r + sphere_c(3);
    sphere_surf = surf(sphere_x, sphere_y, sphere_z);
    set(fig, 'KeyPressFcn', @keypress)

    ii=1;
    for t=0:dt:30
        tic
        
        % Compute shape matching matrices
        A=zeros(d, k, numel(E));
        for i=1:numel(E)
            p = x(:,E{i}) - x0_com;
            Ai = p*B{i};
            A(:,:,i) = Ai;
        end

        % Preparing input for stiffness matrix mex function.
        n=size(x0,2);
        dF_dqij = permute(dF_dq, [3 1 2]);
        Aij = permute(A, [3 1 2]);
        Aij = Aij(:,:);        
        
        % Stiffness matrix (mex function)
        K = -vem3dmesh_neohookean_dq2(Aij, dF_dqij(:,:), ...
                dM_dX(:,:), a, vol, params,k,n,dF,dF_I);
        
        % Force vector
        dV_dq = zeros(numel(x),1);

        % Computing force dV/dq for each point.
        % TODO: move this to C++ :)
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
        
        % Collision force for sphere
        f_collision = zeros(size(x,2)*3,1);
        for i = 1:size(x,2)
          dist_vec = x(:,i) - sphere_c;
          dist = sqrt(sum(dist_vec.^2));
          dir = dist_vec / dist;
          if dist < sphere_r
            f_collision(3*i-2:3*i, 1) = 100 * (sphere_r - dist) * dir;
          end
        end
        f_collision = P * f_collision;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - dt*dt*K)*P') * J;
%         rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error);
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error + f_collision);
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
            fn=sprintf('output/img/cutoff_ten_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end
    
        % Callback to process keypress events
    function keypress(~, evnt)
      switch lower(evnt.Key)  
          case 'leftarrow'
             sphere_c(1) = sphere_c(1) - 10;
          case 'rightarrow'
             sphere_c(1) = sphere_c(1) + 10;
          case 'downarrow'
             sphere_c(2) = sphere_c(2) - 10;
          case 'uparrow'
             sphere_c(2) = sphere_c(2) + 10;
          case 115  % s: move down
             sphere_c(3) = sphere_c(3) - 10;
          case 119  % w: move up
             sphere_c(3) = sphere_c(3) + 10;
          case 97   % a: increase the size of the sphere
             sphere_r = sphere_r - 1;
          case 100  % d: decrease the size of the sphere
             sphere_r = sphere_r + 1; 
          otherwise
              return
      end
      % Always do a redraw
      redraw_sphere();
    end

    function redraw_sphere()
      % draw sphere
      [sphere_x, sphere_y, sphere_z] = sphere;
      sphere_x = sphere_x * sphere_r + sphere_c(1);
      sphere_y = sphere_y * sphere_r + sphere_c(2);
      sphere_z = sphere_z * sphere_r + sphere_c(3);
      sphere_surf.XData = sphere_x;
      sphere_surf.YData = sphere_y;
      sphere_surf.ZData = sphere_z;

      drawnow
    end
end

