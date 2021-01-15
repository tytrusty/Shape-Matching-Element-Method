function quadratic_beam

function pinned_ids = pin_function(x)
    x_max = 1e-2;
    y_min = 0.3;
    y_max = 0.7;
    y_min = 0.4;
    y_max = 0.6;
%     z_min = 0.3;
%     z_max = 0.7;
    z_min = -100.005;
    z_max = 100.5;
    pinned_ids = find(x(1,:) < x_max & x(2,:) > y_min & x(2,:) < y_max ...
                    & x(3,:) > z_min & x(3,:) < z_max);

    pinned_ids1 = find(x(1,:) == min(x(1,:)) & x(2,:) > 0.6 & (x(3,:) > 0.8 | x(3,:) < 0.2));
    pinned_ids2 = find(x(1,:) == min(x(1,:)) & x(2,:) < 0.4 & (x(3,:) > 0.8 | x(3,:) < 0.2));
    pinned_ids = [pinned_ids2 pinned_ids1];
end

iges_file = 'beam2.igs';
% iges_file = 'beam1.igs';
% iges_file = 'beam_1_half.igs';

part = nurbs_from_iges(iges_file);

%material properties
% YM = 5e3; %in Pascals
YM = 5e5; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 2;
options.pin_function = @pin_function;
options.gravity = -50;
options.lambda = lambda;
options.mu = mu;
options.rho = 1;
% options.distance_cutoff = 1.5;
options.distance_cutoff = 0.7;
options.k_stability = 1e5;
% options.k_stability = 1e7;
% options.k_stability = 0;
options.enable_secondary_rays = 0;
options.fitting_mode = 'hierarchical';
options.save_output = 0;

vem_simulate_nurbs(part, options);

end