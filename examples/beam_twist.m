function beam_twist

iges_file = 'beam_cmp3_b.iges';
part = nurbs_from_iges(iges_file);

% cutoff_distance=cutoff_heuristic(part, 0.9,40)

%material properties
YM = 1e6; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt = 0.02;
options.order = 2;
options.pin_function = @pin_function;
options.gravity = 0;
options.lambda = lambda;
options.mu = mu;
options.rho = 1e3;
options.distance_cutoff = 0.51;
options.k_stability = 2e6;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
options.save_obj = 0;
% options.plot_points = 1;
options.x_samples = 15;
options.y_samples = 3;
options.z_samples = 3;
simulate_twist(part, options);
end