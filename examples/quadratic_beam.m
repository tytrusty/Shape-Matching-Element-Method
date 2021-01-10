function quadratic_beam

function pinned_ids = pin_function(x)
    x_max = 1e-2;
    y_min = 0.3;
    y_max = 0.7;
    z_min = 0.3;
    z_max = 0.7;
    z_min = -100.005;
    z_max = 100.5;
    pinned_ids = find(x(1,:) < x_max & x(2,:) > y_min & x(2,:) < y_max ...
                    & x(3,:) > z_min & x(3,:) < z_max);
end

iges_file = 'beam.iges';

part = nurbs_from_iges(iges_file);

options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = 100;
options.mu = 1500;
options.sample_interior = 0; 
options.rho = 1;
options.distance_cutoff = 0.2;

vem_simulate_nurbs(part, options);

end