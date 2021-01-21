function sim=vem_preprocess(parts, varargin)
    % Simulation parameter parsing. Using the same parameters as the
    % simulate method cause I'm lazy and don't want to carry two different
    % options structures.
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
    addParameter(p, 'com_threshold', 100);

    parse(p,varargin{:});
    config = p.Results;
        
    d = 3;              % dimension (2 or 3)
    n = numel(parts);	% number of shapes
    
    % The number of elements in the monomial basis.
    k = basis_size(d, config.order);
    
    % Assembles global generalized coordinates
    [J, hires_J, q, E, x0] = nurbs_assemble_coords(parts);
    
    
    % Sampling points used to compute energies.
    if config.sample_interior
        [V, vol] = raycast_quadrature(parts, [3 3], 10);
    else
    	V=x0;
        vol=ones(size(V,2),1);
    end

    % Compute Shape weights
    [w, w_I] = nurbs_blending_weights(parts, V', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    [w0, w0_I] = nurbs_blending_weights(parts, x0', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    
    % Generate centers of mass.
    [x0_coms, com_cluster, com_map] = generate_com(x0, E, w, n);
    
    % Shape Matrices
    L = compute_shape_matrices(x0, x0_coms, com_map, E, ...
        com_cluster, config.order, config.fitting_mode);

    % Build Monomial bases for all quadrature points
    [Y,Y_S] = vem_dx_dc(V, x0_coms, w, w_I, com_map, config.order, k);
    [Y0,Y0_S] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, config.order, k);
    
    % Computing each gradient of deformation gradient with respect to
    % projection operator (c are polynomial coefficients)
    [dF_dc, dF_dc_S] = vem_dF_dc(V, x0_coms, w, w_I, com_map, config.order, k);

    % Compute mass matrices
    ME = vem_error_matrix(Y0, Y0_S, L, d);
    M = vem_mass_matrix(Y, Y_S, L, config.rho .* vol);
    M = (M + config.k_stability*ME);
    
    % Lame parameters concatenated.
    params = [config.mu * 0.5, config.lambda * 0.5];
    params = repmat(params,m,1);
    
    % Assinging workspace variables to struct.
    sim.k = k;
    sim.J = J;
    sim.hires_J = hires_J;
    sim.q = q;
    sim.E = E;
    sim.x0 = x0;
    sim.V = V;
    sim.vol = vol;
    sim.x0_coms = x0_coms;
    sim.com_cluster = com_cluster;
    sim.com_map = com_map;
    sim.L = L;
    sim.Y = Y; sim.Y_S = Y_S;
    sim.Y0 = Y0; sim.Y0_S = Y0_S;
    sim.dF_dc = dF_dc; sim.dF_dc_S = dF_dc_S;
    sim.w = w; sim.w_I = w_I;
    sim.w0 = w0; sim.w0_I = w0_I;
    sim.ME = ME;
    sim.M = M;
    sim.params = params;
end

