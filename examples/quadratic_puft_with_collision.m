function quadratic_puft_with_collision

function pinned_ids = pin_function(x)
  pinned_ids = [];
end
 
iges_file = 'puft_full_simpler.iges'; % simplified version

sample_density=2;
part = nurbs_from_iges(iges_file,sample_density);

% YM = 2e3; %in Pascals
% pr =  0.25;
% YM = 1e5;
YM = 1e6;
pr = 0.47;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt=0.01;
options.rho = 1.27e3;

options.order = 2;
options.pin_function = @pin_function;
options.gravity = 0;
options.enable_secondary_rays = true;
options.lambda = lambda;
options.mu = mu;
% options.sample_interior = 1;

options.plot_points=0;
options.x_samples=3; %2 on regular 3 on bad
options.y_samples=8;
options.z_samples=15;
% options.z_samples=11;

% New cutoff heuristic thing. The second argument is the "sparsity" which
% is a 0 to 1 value for how local you want deformations to be.
% So 0.9 yields a cutoff distance value that is very local.
cutoff_distance=cutoff_heuristic(part, 0.9);

options.distance_cutoff = cutoff_distance;
% options.distance_cutoff=0.2;
options.save_obj = true;
options.save_obj_path = 'output/puft_with_puft/';
options.k_stability = 1e7;%1e10;

% options.collision_ratio = 4e-1;
options.collision_ratio = 5e4;
options.collision_with_other = true;
options.collision_other_position = [20 0 0];
options.self_collision = false;
options.collision_with_plane = false;
options.collision_plane_z = -20;
options.initial_velocity = [10 0 0];
vem_simulate_nurbs_with_collision(part, options);
  
end