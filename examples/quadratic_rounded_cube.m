function quadratic_rounded_cube

function pinned_ids = pin_function(x)
    % pinning corner points of the cube
    x_min = 0.4;
    y_min = 0.4;
    pinned_ids = find(x(1,:) > x_min & x(2,:) > y_min);
end

iges_file = 'rounded_cube.iges';
part = nurbs_from_iges(iges_file);
YM = 1e6; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);
options.order = 2;
options.rho = 1e3;
options.gravity = -10;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.k_stability=YM;
options.distance_cutoff=cutoff_heuristic(part, 0.9);

% vem_simulate_nurbs_newtons(part, options);
vem_simulate_nurbs(part, options);

end