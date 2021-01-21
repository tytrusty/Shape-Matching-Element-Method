function linear_rocket_with_collision

function pinned_ids = pin_function(x)
    pinned_ids = [];
end

iges_file = 'rocket.iges';

part = nurbs_from_iges(iges_file);

YM = 0.01e9; %in Pascals
pr = 0.5;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 1;
options.gravity = -25.95;
options.rho = 1.3e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 10;

options.k_stability = 1e9;
% 
% options.collision_with_other = false;
% options.self_collision = false;
% options.collision_with_plane = false;
% options.collision_plane_z = -20.0;

vem_simulate_nurbs(part, options);
end