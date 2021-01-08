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
    addParameter(p, 'enable_secondary_rays', true);
    parse(p,varargin{:});
    config = p.Results;
    
    d = 3;  % dimension (2 or 3)
    
    % The number of elements in the monomial basis.
    k = basis_size(d, config.order);
    
    % Read in NURBs 
    fig=figure(1);
    clf;
    
    % resolution = repelem(8,9); resolution(1)=11; % resolution(17)=11;
    parts=nurbs_plot(parts);
    
    % Assembles global generalized coordinates
    [J, ~, q, E, x0] = nurbs_assemble_coords(parts);

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
    [~,L] = compute_shape_matrices(x0, x0_com, E, config.order);
    
    % Build Monomial bases for all quadrature points
    Y = monomial_basis_matrix(V, x0_com, config.order, k);
    Y0 = monomial_basis_matrix(x0, x0_com, config.order, k);
    
    % Compute Shape weights
    w = nurbs_blending_weights(parts, V', config.distance_cutoff, ...
                               config.enable_secondary_rays);
    w_x = nurbs_blending_weights(parts, x0', config.distance_cutoff, ...
                                 config.enable_secondary_rays);
    [W, W_I, W_S] = build_weight_matrix(w, d, k, 'Truncate', true);
    [W0, ~, W0_S] = build_weight_matrix(w_x, d, k, 'Truncate', false);
    
    % Fixed x values.
    x_fixed = zeros(size(x0));
    for i = 1:numel(pin_I)
        x_fixed(:,pin_I(i))=x0(:,pin_I(i));
    end
    
    % Applying fixed point constraints to NURBS jacobian.
    J = P * J;
    
    
    m = size(V,2);  % number of quadrature points
    n = numel(E);	% number of shapes
    
    % Forming gradient of monomial basis w.r.t X
    dM_dX = monomial_basis_grad(V, x0_com, config.order);
    
    % Computing each gradient of deformation gradient with respect to
    % projection operator (c are polynomial coefficients)
    dF_dc = vem_dF_dc(dM_dX, W);

    % Compute mass matrices
    ME = vem_error_matrix(Y0, W0, W0_S, L);
    M = vem_mass_matrix(Y, W, W_S, L);
    M = ((config.rho*M + config.k_stability*ME)); % sparse?
    % Save & load these matrices for large models to save time.
    % save('saveM.mat','M');
    % save('saveME.mat','ME');
    % M = matfile('saveM.mat').M;
    % ME = matfile('saveME.mat').ME;

    ii=1;
    for t=0:config.dt:30
        tic
        
        % Preparing input for stiffness matrix mex function.
        % TODO: don't form this vector this way :)
        b = [];
        for i=1:numel(E)
            b = [b x(:,E{i}) - x0_com];
        end
        b = b(:);

        % Solve for polynomial coefficients (projection operators).
        c = L * b;
        p = c(end-d+1:end);
        x_com = x0_com + p;   
        
        % Stiffness matrix (mex function)
        K = -vem3dmesh_neohookean_dq2(c, dM_dX(:,:), vol, params, ...
                                      dF_dc, W, W_S, W_I, k, n);
        K = L' * K * L;
        
        % Force vector
        dV_dq = zeros(d*(k*numel(E) + 1),1);

        % Computing force dV/dq for each point.
        % TODO: move this to C++ :)
        for i = 1:m
            dMi_dX = squeeze(dM_dX(i,:,:));
            
            % Deformation Gradient
            F = dMi_dX * W{i} * W_S{i} * c;
            F = reshape(F,d,d);

            % Force vector
            dV_dF = neohookean_tet_dF(F, config.lambda, config.mu);

            % Todo: I should be multiplying by volume here, right?
            dV_dq = dV_dq + W_S{i}' * dF_dc{i}' * dV_dF;
        end
        dV_dq = L' * dV_dq;
        
        % Error correction force
        x_centered = x(:);
        x_centered(1:d:end) = x_centered(1:d:end) - x_com(1);
        x_centered(2:d:end) = x_centered(2:d:end) - x_com(2);
        x_centered(3:d:end) = x_centered(3:d:end) - x_com(3);
        f_error = - 2 * ME * x_centered;
        
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
            xi = x(:,x_idx+1:x_idx+x_sz);
            parts{i}.plt.Vertices =xi';
            x_idx = x_idx+x_sz;
        end
        drawnow
        
        if config.save_obj
            obj_fn = "output/obj/part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,parts,obj_fn,ii);
        end
        
        if config.save_output
            fn=sprintf('output/img/cutoff_ten_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end
end

