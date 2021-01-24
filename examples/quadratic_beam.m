function quadratic_beam

function pinned_ids = pin_function(x)
    pinned_ids1 = find(x(1,:) == min(x(1,:)) & x(2,:) > 0.6 & (x(3,:) > 0.8 | x(3,:) < 0.2));
    pinned_ids2 = find(x(1,:) == min(x(1,:)) & x(2,:) < 0.4 & (x(3,:) > 0.8 | x(3,:) < 0.2));
    pinned_ids = [pinned_ids2 pinned_ids1];
end

iges_file = 'beam2.igs';

% To avoid singular nurbs jacobian with excessive pinning, I up the sample
% density to ensure we have enough unpinned samples.
sample_density=3;

part = nurbs_from_iges(iges_file, sample_density);

%material properties
% YM = 5e3; %in Pascals
YM = 5e5; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 2;
options.pin_function = @pin_function;
% options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
options.rho = 1;
% options.distance_cutoff = 1.5;
% options.distance_cutoff = 0.9;
options.distance_cutoff = 1;
options.k_stability = 1e5;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
% options.plot_points = 1;
options.x_samples = 50;
options.y_samples = 1;
options.z_samples = 1;
% vem_simulate_nurbs_newtons(part, options);
vem_simulate_nurbs(part, options);

end