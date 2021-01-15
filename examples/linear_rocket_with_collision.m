function linear_rocket

function pinned_ids = pin_function(x)
    pinned_ids = [];
end

iges_file = 'rocket.iges';

part = nurbs_from_iges(iges_file);

YM = 3e1; %in Pascals
pr =  0.15;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 1;
options.gravity = -25.95;
options.rho = 1e-4;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 10;

% options.collision_with_other = true;
options.self_collision = true;
options.collision_with_plane = true;
options.collision_plane_z = -50.0;
vem_simulate_nurbs_with_collision(part, options);

end