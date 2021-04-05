function octopus

function pinned_ids = pin_function(x)
  pinned_ids = [];
end

% iges_file = 'octopus4.iges';
iges_file = 'squid.iges';

part = nurbs_from_iges(iges_file);

% cutoff_distance=cutoff_heuristic(part, 0.9,40)

% YM = 2e3; %in Pascals
% pr =  0.25;
% YM = 1e5;
% YM = 1.5e4;
YM = 2e4;%2e4;
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
options.save_output=0; % save imgs
options.plot_points=1;
options.plot_com=1;

% Not too bad %
options.x_samples=4;%2; 
options.y_samples=16;
options.z_samples=16;

options.x_samples=6;%4;
options.y_samples=40;
options.z_samples=50;

% options.distance_cutoff=0.70; octupus1
% options.distance_cutoff=0.25;%0.2;% octopus4
options.distance_cutoff=0.08;%0.12;%0.25;%0.2;% octopus4
options.save_obj = 1;
% options.k_stability = YM;%5e3;
options.k_stability = 1e4;% CHANGED 6e3;
% options.k_stability = 2e3;

% changed sampling 8->10, k_stability, and collision ratio 
options.collision_ratio = 0.02;%0.5;
options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = 0.1;
vem_simulate_nurbs_with_collision(part, options);
  
end