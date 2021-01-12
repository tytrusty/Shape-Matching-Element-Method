function linear_trimmed_mug_with_collision

function pinned_ids = pin_function(x)
    verts_to_pin = 3; 
    [~,I] = maxk(x(1,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'mug_complex.iges';

part = nurbs_from_iges(iges_file);

YM = 2e3; %in Pascals
pr =  0.25;
[lambda, mu] = emu_to_lame(YM, pr);

options.rho = 0.1;
options.order = 1;
options.pin_function = @pin_function;
options.gravity = -10;
options.enable_secondary_rays = false;
options.lambda = lambda;
options.mu = mu;
options.sample_interior = 0;
options.collision_ratio = 0.1;
vem_simulate_nurbs_with_collision(part, options);
  
end