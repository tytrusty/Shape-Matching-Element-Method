function linear_mug_implicit

function pinned_ids = pin_function(x)
    pinned_ids = find(x(1,:) > 6 & x(3,:) > 6 & x(3,:) < 7);
end

iges_file = 'mug.iges';

part = nurbs_from_iges(iges_file);

%material properties
YM = 5e2; %in Pascals
pr =  0.15;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 1;
options.gravity = -10;
options.rho = .2;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.sample_interior = 0; % only sample on boundary
options.distance_cutoff = 1;
vem_simulate_nurbs_fmincon_param(part, options);  %%%%% verify hierarchical works!!!!!
                                                        
 
end