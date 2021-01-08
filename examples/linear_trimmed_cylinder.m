function linear_trimmed_cylinder

function pinned_ids = pin_function(x)
    verts_to_pin = 2; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'cylinder.iges';

part = nurbs_from_iges(iges_file);

options.order = 1;
options.pin_function = @pin_function;
options.gravity = -10;
vem_simulate_nurbs(part, options);
  
end