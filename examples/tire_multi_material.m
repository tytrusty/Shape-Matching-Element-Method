function tire_multi_material

function pinned_ids = pin_function(x)
    pinned_ids = [];
end

iges_file = 'tire_5.iges';

part = nurbs_from_iges(iges_file, 1);

YM = 1e6; %in Pascals
% YM = 4e5; %in Pascals
pr = 0.47;
[lambda, mu] = emu_to_lame(YM, pr); 

options.dt = 0.01;
options.order = 2;
options.gravity = -10;
options.rho = 2e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;

% cutoff_distance=cutoff_heuristic(part, 0.9);
options.distance_cutoff = 1.1;%1.2%0.9;%0.7; % 0.75 use 5 for a single com;

options.save_obj = true;

options.k_stability = 1e9; %3e7

% Increase this would give a larger penalty force for collision and
% bounces back to a higher position
options.collision_ratio = 3e2;    % Decrease to 9e5 or 8e5 to make it bounce back lower
options.collision_ratio = 5e3;    % Decrease to 9e5 or 8e5 to make it bounce back lower
% options.collision_ratio = 3e4;    % Decrease to 9e5 or 8e5 to make it bounce back lower
% options.collision_ratio = 4e4;    % Decrease to 9e5 or 8e5 to make it bounce back lower

options.x_samples = 4;
options.y_samples = 9;
options.z_samples = 9; % we want more samples along vertical axis
% If you want centers of mass on the legs of the rocket, you'll need to
% make sure points are sampled on the legs. You can fidget with the above
% numbers to try to change it, if you need.
% options.plot_points=true; % Enable to see all the equadrature points.

options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -7;
options.initial_velocity = [5 0 -5];
options.save_output=0;

vem_tmp_sim(part, options);
end