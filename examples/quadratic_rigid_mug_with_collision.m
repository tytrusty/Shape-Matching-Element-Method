function quadratic_rigid_mug_with_collision

function pinned_ids = pin_function(x)
  pinned_ids = [];
end

% iges_file = 'mug_trimmed.iges';
iges_file = 'mug_trimmed_rot2.iges';

part = nurbs_from_iges(iges_file);

% YM = 2e3; %in Pascals
% pr =  0.25;
% YM = 1e5;
YM = 4e7;
pr = 0.40;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt=0.01;
options.rho = 1e3;

options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.enable_secondary_rays = true;
options.lambda = lambda;
options.mu = mu;
options.sample_interior = 1;

% options.plot_points=1;
options.x_samples=2; %2 on regular 3 on bad
% options.y_samples=5;
% options.y_samples=5;
options.z_samples=15;

options.distance_cutoff=0.1;
options.distance_cutoff=0.3; %0.5;
% options.distance_cutoff=5.6;
options.save_obj = 1;
% options.k_stability = 1e9;

options.initial_velocity = [1.5 0.1 -2];


% options.collision_ratio = 4e-1;
options.collision_ratio = 8e1;
options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -0.8;
vem_simulate_nurbs_with_collision(part, options);
  
end