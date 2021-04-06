function [pcnt,timing]=scalability_test_step(parts,varargin)
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
    addParameter(p, 'save_iges', 0);          	% (0 or 1) whether to output iges files
    addParameter(p, 'save_resultion', 20);    	% the amount of subdivision for the output obj
    addParameter(p, 'pin_function', @(x) 1);
    addParameter(p, 'sample_interior', 1);
    addParameter(p, 'distance_cutoff', 20);
    addParameter(p, 'enable_secondary_rays', true);
    addParameter(p, 'fitting_mode', 'hierarchical');
    addParameter(p, 'plot_points', false);
    addParameter(p, 'plot_com', true);
    addParameter(p, 'initial_velocity', [0 0 0]);
    addParameter(p, 'x_samples', 5);
    addParameter(p, 'y_samples', 9);
    addParameter(p, 'z_samples', 9);
    addParameter(p, 'f_external', [0 0 0]);
    addParameter(p, 'f_external_time', 1000);
    addParameter(p, 'save_obj_path', 'output/obj/');
    
    parse(p,varargin{:});
    config = p.Results;
    
    tic

    d = 3;              % dimension (2 or 3)
    n = numel(parts);	% number of shapes
    k = basis_size(d, config.order);
    
    % Assembles global generalized coordinates
    [J, ~, q, E, x0, S] = nurbs_assemble_coords(parts);
    
    % Initial deformed positions and velocities
    x = x0;
    qdot = reshape(repmat(config.initial_velocity, size(q,1)/3, 1)', [], 1);
    
    pin_I=[];
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    yz_samples = [config.y_samples config.z_samples];
    [V, vol] = raycast_quadrature(parts, yz_samples, config.x_samples);
    m = size(V,2);
    
%     V_plot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',20);
    
    % Lame parameters concatenated.
    params = [config.mu * 0.5, config.lambda * 0.5];
    params = repmat(params,size(V,2),1);
        
    % Compute Shape weights
    [w, w_I] = nurbs_blending_weights(parts, V', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    [w0, w0_I] = nurbs_blending_weights(parts, x0', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    
    % Generate centers of mass.
    [x0_coms, com_cluster, com_map] = generate_com(x0, E, w, n, parts);
%      com_plt = plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
%                         '.','Color','g','MarkerSize',20);
                    
    % Shape Matrices
    L = compute_shape_matrices(x0, x0_coms, com_map, E, ...
        com_cluster, config.order, config.fitting_mode, S);

    % Build Monomial bases for all quadrature points
    [Y,Y_S,C_I] = vem_dx_dc(V, x0_coms, w, w_I, com_map, config.order, k);
    [Y0,Y0_S,C0_I] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, config.order, k);
    
    % Fixed x values.
    x_fixed = zeros(size(x0));
    for i = 1:numel(pin_I)
        x_fixed(:,pin_I(i))=x0(:,pin_I(i));
    end
    
    % Applying fixed point constraints to NURBS jacobian.
    J = P * J;    
    
    % Computing each gradient of deformation gradient with respect to
    % projection operator (c are polynomial coefficients)
    [dF_dc, dF_dc_S] = vem_dF_dc(V, x0_coms, w, w_I, com_map, config.order, k);
    
    % Compute mass matrices
    ME = vem_error_matrix(Y0, L, w0_I, C0_I, d, k, n);
    M = vem_mass_matrix(Y, L, config.rho .* vol, w_I, C_I, d, k, n);
    
    M = (M + config.k_stability*ME); % sparse?
    
    for ii=1:1
        % Preparing input for stiffness matrix mex function.
        b = [];
        for i=1:numel(E)
            b = [b (x(:,E{i}))];
        end
        b = b(:);

        % Solve for polynomial coefficients (projection operators).
        c = L * b;

        % Stiffness matrix (mex function)
        K = -vem3dmesh_neohookean_dq2(c, vol, params, dF_dc, w_I, k, n, ...
                                      size(x0_coms,2));
        K = L' * K * L;
%         figure(3);spy(K);        
%         absK=abs(K);
%         Lz = abs(K(:)) < 1e-8;
%         K(Lz) = 0;
%         [abc,eef]=min(absK(absK>0));
%         fprintf('minval %f',min(absK(absK>0)));
        pcnt=nnz(K)/numel(K);
%         row_nnzs = sum(abs(K) > 1e-8, 2);
%         row_max = max(row_nnzs);
%         row_min = min(row_nnzs);
% 
        % Force vector
        dV_dq = zeros(d*(k*n + size(x0_coms,2)),1);
        
        % Computing force dV/dq for each point.
        % TODO: move this to C++ :)
        for i = 1:m
            % Deformation Gradient
            F = dF_dc{i} * dF_dc_S{i} * c;
            F = reshape(F,d,d);
            
            V(:,i) = Y{i} * Y_S{i} * c;
            
            % Force vector
            dV_dF = neohookean_tet_dF(F, params(i,1), params(i,2));
            dV_dq = dV_dq +  dF_dc_S{i}' * dF_dc{i}' * dV_dF * vol(i);
        end        
        dV_dq = L' * dV_dq;
        
        % Error correction force
        f_error = - 2 * ME * x(:);
        f_error = config.k_stability*(config.dt * P * f_error(:));
       
        % Force from potential energy.
        f_internal = -config.dt*P*dV_dq;
        
        % Computing linearly-implicit velocity update
        % Note: I believe i'm forgetting the error matrix stiffness matrix
        %       but I don't wanna break things so I haven't added it yet.
        lhs = J' * (P*(M - config.dt*config.dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_error);
        qdot = lhs \ rhs;

        % Update position
        q = q + config.dt*qdot;
        x = reshape(P'*J*q,3,[]) + x_fixed;
    end
    timing=toc;

end