function linear_rocket_soft_with_collision

function pinned_ids = pin_function(x)
    pinned_ids = [];
end

iges_file = 'rocket.iges';

part = nurbs_from_iges(iges_file);

YM = 7e6; %in Pascals
pr = 0.4;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 2;
options.gravity = -9.8;
options.rho = 1.27e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 10;

options.save_obj = true;
options.save_obj_path = 'output/obj_rocket_soft_quadratic/';

options.k_stability = YM*1e5;
options.collision_ratio = 1e5;

options.x_samples = 5;
options.y_samples = 9;
options.z_samples = 15; % we want more samples along vertical axis

options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -20.0;
options.initial_velocity = [0 0 -10];

vem_simulate_nurbs_with_collision(part, options);

end