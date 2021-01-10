function linear_trimmed_wrench

function pinned_ids = pin_function(x)
    verts_to_pin = 2; 
    [~,I] = mink(x(1,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'wrench.iges';

part = nurbs_from_iges(iges_file);

options.order = 1;
options.pin_function = @pin_function;
options.gravity = -10;
options.enable_secondary_rays = true;
options.lambda = 0.5 * 1700;
options.mu = 0.5 * 15000;
options.sample_interior = 0;
vem_simulate_nurbs(part, options);
  
end