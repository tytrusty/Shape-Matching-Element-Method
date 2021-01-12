function linear_rocket_with_nose

function pinned_ids = pin_function(x)
    pinned_ids = find(x(3,:) == max(x(3,:)));
end

iges_file = 'rocket_with_nose.iges';
part = nurbs_from_iges(iges_file);

YM = 3e1; %in Pascals
pr =  0.15;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 1;
options.gravity = -25.95;
options.rho = 1e-4;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;

vem_simulate_nurbs(part, options);

end