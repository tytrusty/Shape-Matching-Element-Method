function quadratic_chicken_with_collision

function pinned_ids = pin_function(x)
  pinned_ids = [];
end

% iges_file = 'chicken.iges'; % complex version 
iges_file = 'chicken_simple_closed_small.iges'; % simplified version

part = nurbs_from_iges(iges_file);

% YM = 2e3; %in Pascals
% pr =  0.25;
% YM = 1e5;
YM = 1e4;
pr = 0.47;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt=0.01;
options.rho = 1.27e3;

options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.enable_secondary_rays = true;
options.lambda = lambda;
options.mu = mu;
% options.sample_interior = 1;

options.plot_points=1;
options.x_samples=3; %2 on regular 3 on bad
options.y_samples=6;
options.z_samples=6;
% options.z_samples=11;

% New cutoff heuristic thing. The second argument is the "sparsity" which
% is a 0 to 1 value for how local you want deformations to be.
% So 0.9 yields a cutoff distance value that is very local.
cutoff_distance=cutoff_heuristic(part, 0.9);

options.distance_cutoff = cutoff_distance;
% options.distance_cutoff=0.2;
options.save_obj = true;
options.save_obj_path = 'output/chicken_closed_soft/';
options.k_stability = 1e7;

% options.collision_ratio = 4e-1;
options.collision_ratio = 6e-3;
options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -0.2;
options.initial_velocity = 0.1 * [-10 0 -4];
vem_simulate_nurbs_with_collision(part, options);
  
end