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
options.gravity = -1e3;
options.rho = 8e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 10;

options.save_obj = true;

options.k_stability = 1e9;
options.collision_ratio = 1e6;

options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -20.0;
vem_simulate_nurbs_with_collision(part, options);

end