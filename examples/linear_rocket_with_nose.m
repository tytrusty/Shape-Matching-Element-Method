function linear_rocket_with_nose

function pinned_ids = pin_function(x)
    pinned_ids = find(x(3,:) == max(x(3,:)));
end

iges_file = 'rocket_with_nose.iges';

resolution = repelem(8, 14);

% The first NURBS is the rocket hull, which should receive more samples.
resolution(1) = 13;
resolution(2) = 11;
part = nurbs_from_iges(iges_file, resolution, 1);

options.order = 1;
options.gravity = -1;
options.rho = .1;
options.pin_function = @pin_function;
options.lambda = 170;
options.mu = 1500;

vem_simulate_nurbs(part, options);

end