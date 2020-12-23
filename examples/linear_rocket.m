function linear_rocket

function pinned_ids = pin_function(x)
    pinned_ids = find(x(3,:) < 2);
end

iges_file = 'rocket.iges';

resolution = repelem(7, 13);

% The first NURBS is the rocket hull, which should receive more samples.
resolution(1) = 11;
part = nurbs_from_iges(iges_file, resolution, 1);

options.order = 1;
options.gravity = -20;
options.rho = .2;
options.pin_function = @pin_function;
options.lambda = 100;
options.mu = 1500;

vem_simulate_nurbs(part, options);

end