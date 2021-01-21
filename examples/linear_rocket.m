function linear_rocket

function pinned_ids = pin_function(x)
    pinned_ids = find(x(3,:) < 0.1);
end

iges_file = 'rocket_small.iges';

part = nurbs_from_iges(iges_file);

YM = 3e2; %in Pascals
pr =  0.25;
[lambda, mu] = emu_to_lame(YM, pr);

% options.dt = 0.1;
options.order = 1;
options.gravity = -25.95;
options.rho = 1e-3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 1;
% options.k_stability = YM*1e8;
options.plot_points = 1;
options.x_samples = 2;
options.z_samples = 15;

% vem_simulate_nurbs_newtons(part, options);
vem_simulate_nurbs(part, options);

end