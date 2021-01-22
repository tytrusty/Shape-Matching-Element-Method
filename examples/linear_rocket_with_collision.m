function linear_rocket_with_collision

function pinned_ids = pin_function(x)
    pinned_ids = [];
end

iges_file = 'rocket.iges';

part = nurbs_from_iges(iges_file);

YM = 2e11; %in Pascals
pr = 0.32;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 1;
options.gravity = -20;
options.rho = 8e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;

options.distance_cutoff = 10; % use 5 for a single com;

options.save_obj = false;

options.k_stability = 1e9;
options.collision_ratio = 1e7;

options.x_samples = 5;
options.y_samples = 9;
options.z_samples = 15; % we want more samples along vertical axis

% If you want centers of mass on the legs of the rocket, you'll need to
% make sure points are sampled on the legs. You can fidget with the above
% numbers to try to change it, if you need.
% options.plot_points=true; % Enable to see all the equadrature points.

options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -20.0;
options.initial_velocity = [0 0 -10];

vem_simulate_nurbs_with_collision(part, options);
end