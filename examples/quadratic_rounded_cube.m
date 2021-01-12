function quadratic_rounded_cube

function pinned_ids = pin_function(x)
    % pinning corner points of the cube
    x_min = 40;
    y_min = 55;
    pinned_ids = find(x(1,:) > x_min & x(2,:) > y_min);
end

iges_file = 'rounded_cube.iges';
part = nurbs_from_iges(iges_file);

YM = 25; %in Pascals
pr =  0.15;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 2;
options.rho = 1e-5;
options.gravity = -100;
options.pin_function = @pin_function;
options.gravity = -100;
options.lambda = lambda;
options.mu = mu;

vem_simulate_nurbs(part, options);

end