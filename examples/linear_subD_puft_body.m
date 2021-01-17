function linear_subD_puft
function pinned_ids = pin_function(x)
    verts_to_pin = 2; 
    [~,I] = mink(x(1,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'puft_without_head_subd.iges';
% iges_file = 'puft_subD.iges';

part = nurbs_from_iges(iges_file);

options.order = 1;
options.gravity = -10;
% options.rho = .2;
options.pin_function = @pin_function;
options.lambda = 10*3;
options.mu = 150*3;
% options.sample_interior = 0; % only sample on boundary
% options.distance_cutoff = 1;
% options.enable_secondary_rays = false;
options.save_obj = false;
axis equal
vem_simulate_subd(part, options);

end