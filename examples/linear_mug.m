function linear_mug

function pinned_ids = pin_function(x)
    pinned_ids = find(x(1,:) > 6 & x(3,:) > 6 & x(3,:) < 7);
end

iges_file = 'mug.iges';

resolution = [16 23];
part = nurbs_from_iges(iges_file, resolution, 1);

options.order = 1;
options.gravity = -10;
options.rho = .2;
options.pin_function = @pin_function;
options.lambda = 100;
options.mu = 1500;
options.sample_interior = 0; % only sample on boundary

vem_simulate_nurbs(part, options);

end