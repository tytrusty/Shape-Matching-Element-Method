function vem_simulate_nurbs_fmincon(parts, varargin)
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
    addParameter(p, 'fitting_mode', 'hierarchical');
    addParameter(p, 'plot_points', false);
    addParameter(p, 'plot_com', true);
    addParameter(p, 'x_samples', 5);
    addParameter(p, 'y_samples', 9);
    addParameter(p, 'z_samples', 9);

    parse(p,varargin{:});
    config = p.Results;
    
    d = 3;              % dimension (2 or 3)
    n = numel(parts);	% number of shapes
    
    % The number of elements in the monomial basis.
    k = basis_size(d, config.order);
    
    % Read in NURBs 
    fig=figure(1);
    clf;
    parts=nurbs_plot(parts);

    % Assembles global generalized coordinates
    [J, ~, q, E, x0] = nurbs_assemble_coords(parts);
    
    % Initial deformed positions and velocities
    x = x0;
    qdot=zeros(size(q));
    
    % Setup pinned vertices constraint matrix
    pin_I = config.pin_function(x0);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plotting pinned vertices.
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    % Sampling points used to compute energies.
    if config.sample_interior
        yz_samples = [config.y_samples config.z_samples];
        [V, vol] = raycast_quadrature(parts, yz_samples, config.x_samples);
    else
        V=x0;
        vol=ones(size(V,2),1);
    end
    m = size(V,2);  % number of quadrature points
    
    if config.plot_points
        V_plot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',20);
    end
    
    % Lame parameters concatenated.
    params = [config.mu * 0.5, config.lambda * 0.5];
    params = repmat(params,size(V,2),1);
        
    % Gravity force vector.
  	f_gravity = repmat([0 0 config.rho*config.gravity], size(x0,2),1)';
    f_gravity = config.dt*P*f_gravity(:);

    % Compute Shape weights
    [w, w_I] = nurbs_blending_weights(parts, V', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    [w0, w0_I] = nurbs_blending_weights(parts, x0', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);

    % Generate centers of mass.
    [x0_coms, com_cluster, com_map] = generate_com(x0, E, w, n);
    if config.plot_com
        com_plt = plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
                        '.','Color','g','MarkerSize',20);
        hold on;
    end
    
    % Shape Matrices
    L = compute_shape_matrices(x0, x0_coms, com_map, E, ...
        com_cluster, config.order, config.fitting_mode);

    % Build Monomial bases for all quadrature points
    [Y,Y_S] = vem_dx_dc(V, x0_coms, w, w_I, com_map, config.order, k);
    [Y0,Y0_S] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, config.order, k);

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
    ME = vem_error_matrix(Y0, Y0_S, L, d);
    M = vem_mass_matrix(Y, Y_S, L, config.rho .* vol);
    M = (M + config.k_stability*ME); % sparse?
    % Save & load these matrices for large models to save time.
    % save('saveM.mat','M');
    % save('saveME.mat','ME');
    % M = matfile('saveM.mat').M;
    % ME = matfile('saveME.mat').ME;
    
    PMP = P*M*P';
    
    options = optimoptions('fmincon');
    options.Algorithm = 'trust-region-reflective';
    options.Display = 'none';
    options.SpecifyObjectiveGradient = true;
    options.MaxIterations = 10;
    options.HessianFcn = 'objective';
    %options.CheckGradients = true;
   
    ii=1;
    for t=0:config.dt:30
        tic
        
        % Solve for polynomial coefficients (projection operators).
        c = vem3dmesh_polynomial_coefficients(x, L, E);
        
        % solve for the velocity of the next time step
        qdot = fmincon(@energy, 0*qdot, [], [], [],[], [],[], [], options);
        
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
        
        if config.plot_com
            x_coms = c(d*k*n + 1:end); % extract centers of mass
            com_plt.XData = x_coms(1:d:end);
            com_plt.YData = x_coms(2:d:end);
            com_plt.ZData = x_coms(3:d:end);
        end
        
        if config.plot_points
            V_plot.XData = V(1,:);
            V_plot.YData = V(2,:);
            V_plot.ZData = V(3,:);
        end
        drawnow
        
        if config.save_obj
            obj_fn = "output/obj/part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,parts,obj_fn,ii);
        end
        
        if config.save_output
            fn=sprintf('output/img/simulate_beem_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end  
    
    
%objective needs to return energy, gradient and hessian values
function [e, g, H] = energy(qdot_new)
   % Update position
    q_new = q + config.dt*qdot_new;
    x_new = reshape(P'*J*q_new,3,[]) + x_fixed;
    dt = config.dt;

    % Solve for polynomial coefficients (projection operators).
    c = vem3dmesh_polynomial_coefficients(x_new, L, E);

    neohookean_e =  vem3dmesh_neohookean_q(c, vol, params, dF_dc, dF_dc_S, d);

    e = 0.5*qdot_new'*J'*PMP*J*qdot_new - qdot_new'*J'*PMP*J*qdot + ...
        + config.k_stability * x_new(:)' * ME * x_new(:) + ...
        + 0.5 * neohookean_e + ...
        - qdot_new'*J'*f_gravity;

    if nargout > 1
        g_neohookean = vem3dmesh_neohookean_dq(x_new, c, vol, params, dF_dc, dF_dc_S, ME, L, ...
                                       k, n, d, size(x0_coms,2), config.k_stability);
                    
        g = J' * (PMP*J*(qdot_new - qdot) + ...
            + dt*P*g_neohookean + ...
            - f_gravity);

        if nargout > 2  
            K = vem3dmesh_neohookean_dq2(c, vol, params, dF_dc, w_I, k, n, ...
                                      size(x0_coms,2));
            K = L' * K * L;
          
            H = J' * (PMP + dt*dt*P*K*P') * J;
        end
    end
end

end
