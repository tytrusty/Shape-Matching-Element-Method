function linear_subD_torus
function pinned_ids = pin_function(x)
    verts_to_pin = 2; 
    [~,I] = mink(x(1,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'torus_subd.igs';

part = nurbs_from_iges(iges_file);

options.order = 1;
options.gravity = -10;
% options.rho = .2;
options.pin_function = @pin_function;
% options.lambda = 100;
% options.mu = 1500;
% options.sample_interior = 0; % only sample on boundary
% options.distance_cutoff = 1;
options.enable_secondary_rays = false;
options.save_obj = true;
axis equal
vem_simulate_nurbs(part, options);

end