function quadratic_beam

function pinned_ids = pin_function(x)
%     pinned_ids1 = find(x(1,:) == min(x(1,:)) & x(2,:) > 0.6 & (x(3,:) > 0.8 | x(3,:) < 0.2));
%     pinned_ids2 = find(x(1,:) == min(x(1,:)) & x(2,:) < 0.4 & (x(3,:) > 0.8 | x(3,:) < 0.2));
%     pinned_ids = [pinned_ids2 pinned_ids1];
end

iges_file = 'beam_cmp3.igs';

% To avoid singular nurbs jacobian with excessive pinning, I up the sample
% density to ensure we have enough unpinned samples.
part = nurbs_from_iges(iges_file);

cutoff_distance=cutoff_heuristic(part, 0.9,40)


%material properties
% YM = 5e3; %in Pascals
YM = 5e6; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt = 0.02;
options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
options.rho = 1e3;
% options.distance_cutoff = 1;
% options.distance_cutoff = 0.6;
options.distance_cutoff = 0.6;
% options.k_stability = 1e6;
options.k_stability = 1e6;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
options.save_obj = 1;
% options.plot_points = 1;
options.x_samples = 50;
options.y_samples = 2;
options.z_samples = 3;
% options.x_samples = 25;
% options.y_samples = 1;
% options.z_samples = 1;
% vem_simulate_nurbs_newtons(part, options);
vem_simulate_nurbs(part, options);

end