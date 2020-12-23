function vem_simulate_nurbs(parts, varargin)
    % Simulation parameter parsing
    p = inputParser;
    addParameter(p, 'dt', 0.01);                % timestep
    addParameter(p, 'lambda', 0.5 * 1700);    	% Lame parameter 1
    addParameter(p, 'mu', 0.5 * 15000);        	% Lame parameter 2
    addParameter(p, 'gravity', -300);          	% gravity force (direction is -z direction)
    addParameter(p, 'k_stability', 1e5);       	% stiffness for stability term
    addParameter(p, 'order', 1);              	% (1 or 2) linear or quadratic deformation
    addParameter(p, 'rho', 1);                 	% per point density (currently constant)
    addParameter(p, 'save_output', 0);        	% (0 or 1) whether to output images of simulation
    addParameter(p, 'save_obj', 0);          	% (0 or 1) whether to output obj files
    addParameter(p, 'save_resultion', 20);    	% the amount of subdivision for the output obj
    addParameter(p, 'pin_function', @(x) 1);
    addParameter(p, 'sample_interior', 1);
    addParameter(p, 'distance_cutoff', 20);
    parse(p,varargin{:});
    config = p.Results;

    % Read in NURBs 
    fig=figure(1);
    clf;
    
    % resolution = repelem(8,9); resolution(1)=11; % resolution(17)=11;
    parts=nurbs_plot(parts);
    
    % Assembles global generalized coordinates
    [J, q, E, x0] = nurbs_assemble_coords(parts);

    % Initial deformed positions and velocities
    x = x0;
    qdot=zeros(size(q));
    
    % Setup pinned vertices constraint matrix
    pin_I = config.pin_function(x0);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plot all vertices
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    % Raycasting quadrature as described nowhere yet :)
    if config.sample_interior
        [V, vol] = raycast_quadrature(parts, [9 9], 5);
        vol = vol ./ max(vol); % temporary
    else
    	V=x0;
        vol=ones(size(V,2),1);
    end
    
    % TODO: add option to visualization quadrature points
    % plot3(V(1,:),V(2,:),V(3,:),'.','Color','r','MarkerSize',20);
    
    % Lame parameters concatenated.
    params = [config.lambda, config.mu];
    params = repmat(params,size(V,2),1);
        
    % Gravity force vector.
  	f_gravity = repmat([0 0 config.gravity], size(x0,2),1)';
    f_gravity = config.dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
        
    % Shape Matrices
    %E=cell(1);
    %E{1}=1:size(x0,2);
    [B,~] = compute_shape_matrices(x0, x0_com, E, config.order);
    
    % Build Monomial bases for all quadrature points
    Q = monomial_basis(V, x0_com, config.order);
    Q0 = monomial_basis(x0, x0_com, config.order); 
    
    % Compute Shape weights
    a = nurbs_blending_weights(parts, V', config.distance_cutoff);
    a_x = nurbs_blending_weights(parts, x0', config.distance_cutoff);
    
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
    dM_dX = monomial_basis_grad(V, x0_com, config.order);
    
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
       sum(mm1 < 1e-4);
       dF_I{i} = I';
       SdF{i} = sparse(zeros(numel(I), numel(x)));
       ind=sub2ind(size(SdF{i}), 1:numel(I),I);
       SdF{i}(ind)=1;
       dF{i} = m1;
    end

    % Compute mass matrices
    ME = vem_error_matrix(B, Q0, a_x, d, size(x,2), E);
    M = vem_mass_matrix(B, Q, a, d, size(x,2), E);
    M = ((config.rho*M + config.k_stability*ME)); % sparse?
    % Save & load these matrices for large models to save time.
    % save('saveM.mat','M');
    % save('saveME.mat','ME');
    % M = matfile('saveM.mat').M;
    % ME = matfile('saveME.mat').ME;

    k=4;
    if config.order == 2
        k = 10;
    end

    ii=1;
    for t=0:config.dt:30
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
                dM_dX(:,:), a, vol, params,k, n, dF, dF_I);
        
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
            dV_dF = neohookean_tet_dF(F,config.lambda, config.mu);
            dV_dq = dV_dq + dF_dq(:,:,i)' * dV_dF; % assuming constant area
        end
  
        % Error correction force
        f_error = - 2 * ME * x(:);
        f_error = config.k_stability*(config.dt * P * f_error(:));
       
        % Force from potential energy.
        f_internal = -config.dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - config.dt*config.dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error);
        qdot = lhs \ rhs;

        % Update position
        q = q + config.dt*qdot;
        x = reshape(P'*J*q,3,[]) + x_fixed;
        
        % Update NURBs plots
        x_idx=0;
        for i=1:numel(parts)
            x_sz = size(parts{i}.x0,2);
            xi = reshape(x(:,x_idx+1:x_idx+x_sz), 3, parts{i}.subd(1), parts{i}.subd(2));
            parts{i}.plt.XData = squeeze(xi(1,:,:));
            parts{i}.plt.YData = squeeze(xi(2,:,:));
            parts{i}.plt.ZData = squeeze(xi(3,:,:));
            x_idx = x_idx+x_sz;   
        end
        drawnow
        
        if config.save_obj
            obj_fn = "output/obj/part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,parts,config.obj_res,obj_fn,ii);
        end
        
        if config.save_output
            fn=sprintf('output/img/cutoff_ten_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end
end

