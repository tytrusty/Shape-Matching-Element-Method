function linear_trimmed_mug_with_collision

function pinned_ids = pin_function(x)
% %     verts_to_pin = 3; 
%     verts_to_pin = 10; 
%     [~,I] = maxk(x(1,:),verts_to_pin);
%     pinned_ids = I(1:verts_to_pin);
  pinned_ids = [];
end

iges_file = 'mug_complex.iges';

part = nurbs_from_iges(iges_file);

% YM = 2e3; %in Pascals
% pr =  0.25;
YM = 2.9e7;
pr = 0.32;
[lambda, mu] = emu_to_lame(YM, pr);

options.rho = 0.1;
options.order = 1;
options.pin_function = @pin_function;
options.gravity = -100;
options.enable_secondary_rays = false;
options.lambda = lambda;
options.mu = mu;
options.sample_interior = 0;

% options.save_obj = true;

options.collision_ratio = 0.1;
options.collision_with_other = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = -40.0;
vem_simulate_nurbs_with_collision(part, options);
  
end