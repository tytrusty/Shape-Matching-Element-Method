function quadratic_rounded_cube

function pinned_ids = pin_function(x)
    % pinning corner points of the cube
    x_min = 40;
    y_min = 55;
    pinned_ids = find(x(1,:) > x_min & x(2,:) > y_min);
end

iges_file = 'rounded_cube.iges';

% Resolution indicates how many point samples we will take on each
% e.g. 6 means we have 6 samples in both the U & V coordinates, so
%      a total of 36 samples across the NURBs patch.
resolution = 6;
part = nurbs_from_iges(iges_file, resolution, 0);

options.order = 2;
options.pin_function = @pin_function;
options.gravity = -100;
options.lambda = 1700;
options.mu = 15000;
options.rho = 1;

vem_simulate_nurbs(part, options);

end